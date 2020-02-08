import java.sql.*;

/** this application uses the dumpCluserPrt Object to 
	dump the title.
 */

public class dumpClusterTitle {

	public static void main(String[] args) {
		//String pgQuery =  "select pcluster, trim(division) as div, prt from cds join organism on source=org_id where pcluster notnull order by pcluster, division LIMIT 117000 offset 74587";
		String pgQuery =  "select pcluster, trim(division) as div, prt from cds join organism on source=org_id where pcluster notnull order by pcluster, division";

		boolean doThread = true;
		int size = 2000;
		int i = 0;
		while (i < args.length) {
			if (args[i].equals("-T"))  doThread = false; 
			else if (args[i].equals("-s")) size = Integer.parseInt(args[++i]);
			i++;
		}

		try {
			dumpClusterPrt dumper = new dumpClusterPrt(pgQuery);
			if (!doThread) dumper.dumpClusterDef(size);
			else dumper.threadedTitleDump(size);
		}
		catch(ClassNotFoundException ce) {
			System.out.println(ce);
			ce.printStackTrace();
		}
		catch(SQLException se) {
			System.out.println(se);
			se.printStackTrace();
		}
	}

}
