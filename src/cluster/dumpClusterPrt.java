import java.util.*;
import java.sql.*;
import java.io.*;
import acedb.*;
import fetchace.*;

/** This program will dump all the proteins within one cluster into
 * each file.  All the clusters from orpara database stored in the
 * cds table will be dumped.
 */

public class dumpClusterPrt extends Thread
{
	private static String user = "kzhou";
	private static String passwd = "fugufish";
	private static String pgPasswd = ""; // no thing needed
	private static String aceHost = "localhost";
	private static int[] ports = {3001, 3002, 3003, 3004, 3005, 3006, 3007, 3008, 3009};
	/** inside the postgresdb division is char(30), so they must be trimmed in
	 * the query, otherwise the HashMap div2db will not find it!. */
	private static String[] dbs = {"amp", "fr", "hs", "mam", "mm", "pri", "rod", "sau", "vrt"};
	public static File dir = null;
	//private static Ace acedb[8];
	private HashMap div2db = new HashMap();

	private static String pgdb = "jdbc:postgresql://localhost/parortho";
	ResultSet clusters = null;



	/** initialize the 9 Acedb on localhost.  They could be on different machines
	 * for parallel processing */
	public dumpClusterPrt(String pgQuery) throws ClassNotFoundException, SQLException {
		for (int i=0; i<9; i++) {
			//System.out.println(dbs[i] + " " + ports[i]); //debug
			div2db.put(dbs[i], new Ace(aceHost, ports[i], user, passwd));
		}
		Class.forName("org.postgresql.Driver");
		Connection pgconn = DriverManager.getConnection(pgdb, user, pgPasswd);
		Statement stmt = pgconn.createStatement();
		clusters = stmt.executeQuery(pgQuery);
	}

	/** obtains the fastaFormated protein sequence.*/
	public String getPrt(String div, String prtName) throws AceException {
		return ((Ace)div2db.get(div)).getFastaPeptide(prtName);
	}

	/** Get the title or definition for protein prtName. 
	 * @return null if prrName has no Title. */
	public String getDef(String div, String prtName) throws AceException {
		Aceobj aceprt = ((Ace)div2db.get(div)).fetch("Protein", prtName);
		if (aceprt == null) {
			System.err.println("not getting this object: " + prtName + " from " + div);
			System.exit(1);
		}
		aceprt = aceprt.at("Title");
		if (aceprt != null) {
			return aceprt.right().toString();
		}
		else return null;
	}

	/** dump all clusters into files with filesize clusters in each file.
	 * Files will be named 1.def 2.def ... */
	public void dumpClusterDef(int filesize) {
		try {
			String clusterid, div, prt, pep;
			boolean more;
			int clusterCnt=0, dirCnt;
			//PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter("cluster.def")));
			PrintWriter out = null;
			File deffile = null;

			if (!clusters.next()) { System.exit(1); } 
			clusterid = clusters.getString(1);
			div = clusters.getString(2);
			prt = clusters.getString(3);
			while (true) {
				if (clusterCnt%filesize == 0) {
					int fc = clusterCnt/filesize + 1;
					String fn = Integer.toString(fc) + ".def";
					deffile = new File(fn);
					deffile.createNewFile();
					out = new PrintWriter(new BufferedWriter(new FileWriter(deffile)));
				}
				clusterCnt++;

				String title = getDef(div, prt);
				out.println("Cluster: " + clusterid);
				if (title != null) out.println(prt + "->" + title);
				else out.println(prt + "->not found");
				System.out.println("fetching cluster: " + clusterid + " ...");
				while ((more = clusters.next()) && clusterid.equals(clusters.getString(1))) {
					//System.out.println(title);  //debug
					prt = clusters.getString(3);
					title = getDef(clusters.getString(2), prt);
					if (title != null) out.println(prt + "->" + title);
					else out.println(prt + "->not found");
				}
				out.println("==========\n");
				if (!more) break;
				clusterid = clusters.getString(1);
				div = clusters.getString(2);
				prt = clusters.getString(3);
			}
			out.close();
		}
		catch(Exception exc) {
			System.err.println("Exception\n" + exc);
			exc.printStackTrace();
		}
	}

	/** dump the title of each protein sequence into files with 
	 * filesize clusters in each file.  Default filesize=2000. */
	public void threadedTitleDump(int filesize) {
		try {
			String clusterid, div, prt, pep;
			boolean more;
			int clusterCnt=0, dirCnt;
			PrintWriter out = null; 
			File deffile = null;

			if (!(more=clusters.next())) { System.exit(1); } 
			clusterid = clusters.getString(1);
			div = clusters.getString(2);
			while (true) {
				///////////////////////////// making files /////////////////////
				if (clusterCnt%filesize == 0) {
					int fc = clusterCnt/filesize + 1;
					String fn = Integer.toString(fc) + ".def";
					deffile = new File(fn);
					deffile.createNewFile();
					if (out != null) out.close();
					out = new PrintWriter(new BufferedWriter(new FileWriter(deffile)));
				}
				clusterCnt++;
				////////////////////////////////////////////////////////////////////
				out.println("\ncluster: " + clusterid);
				System.err.println("\nWorking on cluster: " + clusterid + " ...");
				//int requestCount = 0;
				while (more && clusterid.equals(clusters.getString(1))) {  // each cluster
					Vector prtnames = new Vector();
					while (more && clusterid.equals(clusters.getString(1)) &&
							               div.equals(clusters.getString(2)) ) { // each division
						prtnames.add(clusters.getString(3));
						//requestCount++;
						more = clusters.next();
					}
					new Fetchacedef((Ace)div2db.get(div), div, prtnames, out).start();
					if (!more) break;
					div = clusters.getString(2);
				}

				while (Fetchacedef.threadCount > 0) {
					//sleep(requestCount);
					sleep(5);
				}
				if (!more) break;
				clusterid = clusters.getString(1);
			}
			out.close();
		}
		catch(Exception exc) {
			System.err.println("Exception\n" + exc);
			exc.printStackTrace();
		}
	}

	/** a helper function to manage files and directories. */
	static PrintWriter createPR(int cnt, String ext) throws IOException {
		int dirCnt;
		if (cnt%2000 == 0) {
			dirCnt = cnt/2000 + 1;
			dir= new File("dir" + Integer.toString(dirCnt));
			dir.mkdir();
		}
		String fileName = Integer.toString(cnt + 1) + "." + ext;
		File outfile = new File(dir, fileName);
		outfile.createNewFile();
		return new PrintWriter(new BufferedWriter(new FileWriter(outfile)));
	}

	public static void main(String[] args) {

		/////////// postgres database connection information ///////////////
		//
		String pgPasswd = "";
		/** the trim(division) is essential to remove to whitespaces */
		String pgQuery = "select pcluster, trim(division) as div, prt from cds join organism on source=org_id where pcluster notnull order by pcluster, division";

		if (args.length == 0) {
			System.out.println(" Usage: java dumpClusterPrt -u user\n -w passowrd\n" );
			System.out.println("default user=" + user);
			//System.exit(1);
		}
		else {
			int i = 0;
			while (i < args.length) {
				if (args[i].equals("-u")) user = new String(args[++i]);
				else if (args[i].equals("-w")) passwd = new String(args[++i]);
				else {
					System.out.println(args[i] + " is an unknow argument");
					System.exit(1);
				}
				i++;
			}
		}


		try {
			dumpClusterPrt dumper = new dumpClusterPrt(pgQuery);

			/*
			Class.forName("org.postgresql.Driver");
			Connection pgconn = DriverManager.getConnection(pgdb, user, pgPasswd);
			Statement stmt = pgconn.createStatement();
			ResultSet rs = stmt.executeQuery(pgQuery);
			*/

			String clusterid, div, prt, pep;
			boolean more;
			int fileCnt=0, dirCnt;
			//File dir = null, outfile = null;
			PrintWriter out;
			//String fileName;
			//File dir = null;

			if (!dumper.clusters.next()) { System.exit(1); } 
			clusterid = dumper.clusters.getString(1);
			div = dumper.clusters.getString(2);
			prt = dumper.clusters.getString(3);
			while (true) {
				out = createPR(fileCnt, "pep");

				/*
				if (fileCnt%2000 == 0) {
					dirCnt = fileCnt/2000 + 1;
					dir= new File("dir" + Integer.toString(dirCnt));
					dir.mkdir();
				}
				fileName = Integer.toString(fileCnt + 1) + ".pep";
				outfile = new File(dir, fileName);
				outfile.createNewFile();
				out = new PrintWriter(new BufferedWriter(new FileWriter(outfile)));
				*/

				fileCnt++;
				pep = dumper.getPrt(div, prt);
				out.print(pep);
				System.out.println("fetching cluster: " + clusterid + " ...");
				while ((more = dumper.clusters.next()) && clusterid.equals(dumper.clusters.getString(1))) {
					pep = dumper.getPrt(dumper.clusters.getString(2), dumper.clusters.getString(3));
					out.print(pep);
				}
				out.close();
				if (!more) break;
				clusterid = dumper.clusters.getString(1);
				div = dumper.clusters.getString(2);
				prt = dumper.clusters.getString(3);
			}
		}
		catch(ClassNotFoundException ce) {
			System.out.println("Database driver problem");
		}
		catch(SQLException se) {
			System.out.println("Querying parortho database failed");
		}
		catch(IOException ie) {
			System.err.println("File openning problem" + ie);
			System.exit(1);
		}
		catch(AceException acex) {
			System.err.println("Ace exec error" + acex);
			System.exit(1);
		}
		catch(Exception exc) {
			System.err.println("Exception\n" + exc);
			exc.printStackTrace();
		}
	}
}
