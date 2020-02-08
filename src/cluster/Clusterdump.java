import java.util.*;
import java.sql.*;
import java.io.*;
import acedb.*;
//import fetchace.Fetchaceprt;
//This has moved to be in the same compilation unit

/** This program will dump all the proteins within one cluster into
 * each file.  All the clusters from orpara database stored in the
 * cds table will be dumped.
 * We put all the sequences from each cluster into one file. We only put certain
 * number of files in each directory, so that we don't get too many files in 
 * each directory.
 *
 *	We put all the title into one file.  Title of each cluster is put next to
 *	each other separateted by empty line.
 */

public class Clusterdump extends Thread
{
	private DBinfo dbinfo = new DBinfo();

	/* This map holds all the div-to->acedb mapping.
	 */
	private HashMap div2db = null; //dbinfo.createAllAce();
	ResultSet result = null;     // holds the query result

	private Ace getAce(String dbname) {
		return (Ace)div2db.get(dbname);
	}

	/** Construct a dump object that will dump the pgQuery result into files of
		protein sequences in fasta format and file of protein titles for further
		analysis.  

		This object may perform threaded dump or simple loop through dump depending
		on the dump method chosen.
	 */
	public Clusterdump(String pgQuery) throws ClassNotFoundException, SQLException {
		div2db = dbinfo.createAllAce();
		Connection pgconn = dbinfo.createPg();
		Statement stmt = pgconn.createStatement();
		result = stmt.executeQuery(pgQuery);
	}

	/** Simple loop through dump both protein fasta sequence and title will
		be dumpped.
		The format of the query is fixed:
		    select pcluster, trim(db) as div, prt 
			 from cds c join species s on c.source=s.id 
			 where pcluster notnull 
			 order by pcluster, db, prt
	  Files will be named 1.def 2.def ... 
	  numf numfile per directory
	 */
	public void dump(int numf, String prefix) {
		try {
			String clusterid, div, prt, pep, title;
			boolean more = true;
			int clusterCnt=0, dirCnt;
			PrintWriter seqout = null;  // sequence file, one file per cluster
			PrintWriter defout = null;  // definition file, one file per directory

			more = result.next();  // start reading the first tuple of the result
			if (!more) { 
				System.err.println("query result is empty");
				System.exit(1); 
			} 
			File subdir = null;

			// loops one cluster at a time
			while (more) {
				/// controls directory creation, and new deffile creation/////////
				if (clusterCnt%numf == 0) {
					dirCnt = clusterCnt/numf + 1;  // count from 1
					String subdirstr = Integer.toString(dirCnt);
					subdir = new File("clusterpepdir" + subdirstr);
					subdir.mkdir();
					File deffile = new File("clusterdef" + subdirstr);  // def file in current dir
					deffile.createNewFile();

					if (defout != null) defout.close();  // not the first time
					defout = new PrintWriter(new BufferedWriter(new FileWriter(deffile)));
				}
				// start reading the tuples from the database
				clusterid = result.getString(1);
				clusterCnt++;

				defout.println("Cluster: " + clusterid);
				System.out.println("fetching cluster: " + clusterid + " ...");
				//while ((more=result.next()) && clusterid.equals(result.getString(1))) {

				//// construct output stream for peptide sequences /////////
				File clusterfile = new File(subdir, prefix + clusterid + ".pep");
				seqout = new PrintWriter(new BufferedWriter(new FileWriter(clusterfile)));
				// loop through all db in this cluster
				while (more && clusterid.equals(result.getString(1))) { // all db of cluster
					div = result.getString(2);
					Ace ace = getAce(div);
					System.err.print(".");  // progres report
					// loop through all members in each db
					while (more && div.equals(result.getString(2))) {
						prt = result.getString(3);
						// System.err.println("Fetching protein: " + prt + "...");

						///////// get title /////////////////
						Aceobj aceprt = ace.fetch("Protein", prt);
						if (aceprt == null) {
							System.err.println("null object: " + prt + " from " + div);
							System.exit(1);
						}
						aceprt = aceprt.at("Title");
						if (aceprt == null) {
							title = "No Title";
						}
						else title = aceprt.right().toString();
						///////// get pep /////////////////////
						pep = ace.getFastaPeptide(prt);
						// System.err.println(pep); // debug

						/// writing to output stream /////////
						defout.println(prt + "->" + title);
						seqout.print(pep);  // no title line for fasta seq
						more = result.next();  // reading the next tuple
					}
				}
				defout.println("==========\n");  // one cluster is done
				seqout.close();
			}
			defout.close();
		}
		catch(Exception exc) {
			System.err.println("Exception\n" + exc);
			exc.printStackTrace();
		}
	}

	/** Exact copy of dump, with thread enabled.  May speed up by 2-100 fold
	 * depends on the setting of servers and clients and the comunication
	 *	 speed of the connection.
    *
    * Result dump into directories and files in the current directory
    * with predefined names.
	 * @param numf number of sequence files per directory
	 * @param prefix is a prefix to give the cluster names, default C.  This is essential for clustalw program.  It got a bug when feeding it with 1.pep file it fails.  with this set, 1.pep will becom C1.pep
	 */
	public void tdump(int numf, String prefix) {
		try {
			String clusterid, div, prt, pep, title;
			boolean more = true;
			int clusterCnt=0, dirCnt;
			PrintWriter seqout = null;  // sequence file, one file per cluster
			PrintWriter defout = null;  // definition file, one file per directory
			more = result.next();  // start reading the first tuple of the result
			File subdir = null;
			
			while (more) { // loops one cluster at a time
				/// controls directory creation, and new deffile creation/////////
				if (clusterCnt%numf == 0) {
					dirCnt = clusterCnt/numf + 1;  // count from 1
					String subdirstr = Integer.toString(dirCnt);
					subdir = new File("clusterpepdir" + subdirstr);
					subdir.mkdir();
					File deffile = new File("clusterdef" + subdirstr);  // def file in current dir
					deffile.createNewFile();

					if (defout != null) defout.close();  // not the first time
					defout = new PrintWriter(new BufferedWriter(new FileWriter(deffile)));
				}
				clusterid = result.getString(1);
				clusterCnt++;

				defout.println("Cluster: " + clusterid);
				System.out.println("fetching cluster: " + clusterid + " ...");

				//// construct output stream for peptide sequences /////////
				File clusterfile = new File(subdir, prefix + clusterid + ".pep");
				seqout = new PrintWriter(new BufferedWriter(new FileWriter(clusterfile)));
				// loop through all db in this cluster
				while (more && clusterid.equals(result.getString(1))) { // all db of cluster
					div = result.getString(2);
					//Ace ace = getAce(div);
					//System.err.print(".");  // progres report
					// loop through all members in each db
					Vector prtlist = new Vector();
					while (more && clusterid.equals(result.getString(1)) &&
							                div.equals(result.getString(2))) 
					{
						prt = result.getString(3);
						prtlist.add(prt);
						more = result.next();  // reading the next tuple
					}
					Fetchaceprt fetcher = new Fetchaceprt(getAce(div), div, prtlist, seqout, defout);
					fetcher.start();
				}
				long wait = 1;  // main thread sleep time
				while (Fetchaceprt.threadCount > 0) {
					Thread.sleep(wait++);
					// Linear incremental waiting policy
					//System.err.println("wait: " + wait); // debug
				}
				defout.println("==========\n");  // one cluster is done
				seqout.close();
			}
			defout.close();
		}
		catch(Exception exc) {
			System.err.println("Exception\n" + exc);
			exc.printStackTrace();
		}
	}

	/** dumps only one cluster given p_cluster = clid.
	 * In this case, the pgQuery string should be 
	 * 	Select db, prt from cds join species s on cds.source=s.id
	 * 	where p_cluster=clusterid order by db, prt
	 *  Where should the result go?
	 *  To files first
	 */
	public void tdumpone(String clusterid) throws IOException, SQLException, InterruptedException {
		String db, prt, pep, title;
		PrintWriter seqout = null;  // sequence file, one file per cluster
		PrintWriter defout = null;  // definition file, one file per directory

		File deffile = new File("cluster" + clusterid + ".def");  // def file in current dir
		deffile.createNewFile();
		defout = new PrintWriter(new BufferedWriter(new FileWriter(deffile)));

		defout.println("Cluster: " + clusterid);

		//// construct output stream for peptide sequences /////////
		File clusterfile = new File("cluster" + clusterid + ".pep");
		seqout = new PrintWriter(new BufferedWriter(new FileWriter(clusterfile)));

		boolean more = result.next();  // start reading the first tuple of the result
		// loop through all db in this cluster
		while (more) { // all db of cluster
			db = result.getString(1);
			Vector prtlist = new Vector();
			// loop through all members in each db
			while (more && db.equals(result.getString(1))) {
				prt = result.getString(2);
				prtlist.add(prt);
				more = result.next();  // reading the next tuple
			}
			Fetchaceprt fetcher = new Fetchaceprt(getAce(db), db, prtlist, seqout, defout);
			fetcher.start();
		}
		//long wait = 1;  // main thread sleep time
		while (Fetchaceprt.threadCount > 0) {
			//Thread.sleep(wait++);  // this takes longer, because this is not a 
			// CPU intensive process, we got extra CPU
			// Linear incremental waiting policy
			Thread.sleep(1);  // fixed waiting shortest waiting time, fastest
			//System.err.println("wait: " + wait); // debug
		}
		seqout.close();
		defout.close();
	}

	public static void main(String[] args) {
	/** should be able to dump either the whole cluster database or
		just dump one cluster.  Right now I have only implemented dumping
		the whole database.
	 */

		/** the trim(division) is essential to remove to whitespaces */
		String pgQuery = "select p_cluster, trim(db) as div, prt from cds c join species s on c.source=s.id where p_cluster notnull order by p_cluster, db, prt";

		boolean thread = true;
		String cluster = null;
		int i=0;
		while (i<args.length) {
			if (args[i].equals("-T"))  thread=false;  // non_thread dump
			else  cluster = args[i];

			/*
			else {
				System.err.println(args[i] + " not a legal argument");
				System.err.println("Usage: java Clusterdump -c cluster_id -T [no thread]");
				System.exit(1);
			}
			*/
			i++;
		}

		/** for one cluster fetch */
		String pgQueryone = "select trim(db) as db, prt from cds c join species s on c.source=s.id where p_cluster=" + cluster + " order by db, prt";
	
		try {
			java.util.Date d1 = new java.util.Date();  // starting time
			Clusterdump dumper;
			if (cluster == null) {
				dumper = new Clusterdump(pgQuery);
				if (thread) dumper.tdump(2000, "C");  // adding prefix C1.pep
				else dumper.dump(2000, "C");
			}
			else {
				System.err.println("fetch all proteins from cluster " + cluster);
				dumper = new Clusterdump(pgQueryone);
				dumper.tdumpone(cluster);  // only the treaded version
			}
			java.util.Date d2 = new java.util.Date(); // ending time
			System.err.println("start: " + d1 + " end: " + d2);
		}
		catch(ClassNotFoundException e1) {
			System.err.println("Class not found for postgres driver");
			e1.printStackTrace();
		}
		catch (SQLException e2) {
			System.err.println("SQL exception");
			e2.printStackTrace();
		}
		catch (InterruptedException e3) {
			System.err.println("thread interrupted\n");
			e3.printStackTrace();
		}
		catch (IOException e4) {
			System.err.println("IO problem");
			e4.printStackTrace();
		}
	}
}

/** fetcheas both the protein part and the title part simultaneously.
	This class should be used in place of Fetchacedef.
 */
class Fetchaceprt extends Thread {
	private Ace acedb = null;
	private Vector prtlist = new Vector();
	private String divname = null;
	private PrintWriter pepout, defout;

	public static int threadCount = 0;  //to keep track of the works done on different thread
	public static final int LARGE = 10; // cluster with more than LARGE members will
	// be fetched with the bulk method.

	/** Constructor.  div is redundant with a.  list is the input 
	 * protein ids. There two outputs pw_pep for protein, pw_def 
	 * for definition.
	 */
	public Fetchaceprt(Ace a, String div, Vector list, PrintWriter pw_pep,
			                                                PrintWriter pw_def) { 
		acedb = a;
		prtlist = list;
		divname = div;
		pepout = pw_pep;
		defout = pw_def;
		increaseCount();
	}

	private synchronized void increaseCount() { ++threadCount; }
	private synchronized void decreaseCount() { --threadCount; }

	/** not used in run but can be used for other purposes. */
	private synchronized String getDef(String prtname)
				     throws AceException {
		Aceobj ptr = acedb.fetch("Protein", prtname);
		if (ptr == null) {
			System.err.println(prtname + " not found in this db");
			System.exit(1);
		}
		ptr = ptr.at("Title");
		if (ptr == null) return null;
		return ptr.right().toString();
	}
	/** not used in run but can be used in other methods.*/
	private synchronized String getPep(String name) throws AceException,
		AceObjectNotFoundException {
		String pep = acedb.getFastaPeptide(name);
		System.out.println("getting protein: " + name);  // debug
		//return acedb.getFastaPeptide(name);
		System.out.println("the prt seq is " + pep);
		return pep;
	}

	/** simple dump not to be used under threaded condition.
	 * this method is the same except that it is not executed under the 
	 * thread environment.
	 */
	public void dump() throws AceException {
		String name = null;
		String title = null;
		for (int i=0; i<prtlist.size(); i++) {
			name = (String)prtlist.get(i);
			Aceobj prt = acedb.fetch("Protein", name);
			if (prt == null) {
				System.err.println(name + " not found in this db");
				System.exit(1);
			}
			prt = prt.at("Title");
			if (prt == null)  title = "not found"; 
			else title = prt.right().toString();
			defout.println(divname + " | " + name + ": " + title);

			String prtseq = acedb.getFastaPeptide(name);
			pepout.print(prtseq);  // seq follow each other without space
		}
	}
	/** create several keyset commands if prtlist is too long with 
	 * each String shorter than Ace.QUERY_LIMIT.
	 */
	public static Vector prepareKeyset(Vector prtlist) {
		Vector trunk = new Vector(10);
		int i=-1;
		while (i+1 < prtlist.size()) {
			StringBuffer tmp = new StringBuffer(Ace.QUERY_LIMIT);
			tmp.append("Keyset-read = Protein ");
			while (i+1 < prtlist.size() && tmp.length()+((String)prtlist.get(i+1)).length() < Ace.QUERY_LIMIT) {
				tmp.append((String)prtlist.get(i+1) + ";");
				i++;
			}
			trunk.add( new String(tmp) );
		}
		return trunk;
	}

	/** Does both title and sequence fetch and dump. */
	public void run() {
		//String name = null;
		//String title = null;
		Aceobj t = null;
		try {
			if (prtlist.size() > LARGE) {
				Vector keysetTrunks = prepareKeyset(prtlist);
				for (int i=0; i<keysetTrunks.size(); i++) {
					acedb.keyset((String)keysetTrunks.get(i));
					String allpep = acedb.exec("Peptide");
					pepout.print(Ace.chopComment(allpep));
					//System.out.println(allpep);  // debug
					//System.out.println(acedb.fetchManyString("Protein", prtlist)); //debug
					Objset allprt = acedb.fetchMany((String)keysetTrunks.get(i));
					Aceobj prt_obj = null;
					while ((prt_obj = allprt.next()) != null) {
						t = prt_obj.at("Title");
						if (t != null) t = t.right();
						defout.println(divname + " | " + prt_obj + ": " + t);
					}
				}
			}
			else {
				for (int i=0; i<prtlist.size(); i++) {
					Aceobj prt = acedb.fetch("Protein", (String)prtlist.get(i));
					if (prt == null) {
						System.err.println(prtlist.get(i) + " not found in this db");
						System.exit(1);
					}
					t = prt.at("Title");
					if (t != null)  t = t.right(); 
					defout.println(divname + " | " + prt + ": " + t);

					pepout.print( acedb.getFastaPeptide((String)prtlist.get(i)) );
				}
			}
			decreaseCount();
		}
		catch (AceException ae) {
			System.err.println("!!!! Ace problem " + ae);
			ae.printStackTrace();
			threadCount--;
			System.exit(1);
		}
		catch (Exception e) {
			System.err.println("!!!! Fetal problem " + e);
			threadCount--;
			e.printStackTrace();
			System.exit(1);
		}
	}
}
