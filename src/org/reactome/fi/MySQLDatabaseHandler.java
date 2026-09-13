/*
 * Created on Jul 10, 2012
 *
 */
package org.reactome.fi;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.Statement;
import java.util.HashMap;
import java.util.Map;

import org.apache.log4j.Logger;
import org.gk.persistence.MySQLAdaptor;
import org.junit.Test;
import org.reactome.fi.util.FIConfiguration;

/**
 * This class is used to dump mysql database and convert the dumped file from Innodb to MyIsam.
 * The conversion from Innodb to MyISAM cannot be done for the FI database because of foreign keys
 * constraints imposed by hibernate.
 * @author gwu
 *
 */
public class MySQLDatabaseHandler {
    private static final Logger logger = Logger.getLogger(MySQLDatabaseHandler.class);
    
    public MySQLDatabaseHandler() {
    }
    
    
    /**
     * Under Mac, the table names are case insensitive, but under Linux, the table names are case sensitive.
     * In some occasions, the table names under Mac are in lower case. This method is used to fix the table names.
     * @throws Exception
     */    
    @Test
    public void fixTableNameCases() throws Exception {
        // This is the previous version, configured via PREV_REACTOME_SOURCE_DB_NAME
        // (derived from PREV_REACTOME_RELEASE in build.params) instead of a hardcoded
        // release number, since that changes every year.
        MySQLAdaptor adapter = new MySQLAdaptor("localhost",
                FIConfiguration.getConfiguration().get("PREV_REACTOME_SOURCE_DB_NAME"),
                FIConfiguration.getConfiguration().get("DB_USER"),
                FIConfiguration.getConfiguration().get("DB_PWD")
                );
        Map<String, String> tableNameMap = new HashMap<>();
        Connection conn = adapter.getConnection();
        Statement stmt = conn.createStatement();
        ResultSet rs = stmt.executeQuery("SHOW TABLES");
        while (rs.next()) {
            String tableName = rs.getString(1);
            String lowerCaseTableName = tableName.toLowerCase();
            tableNameMap.put(lowerCaseTableName, tableName);
        }
        rs.close();
        stmt.close();
        System.out.println("Table names in the previous version ("
                + FIConfiguration.getConfiguration().get("PREV_REACTOME_SOURCE_DB_NAME") + "): "
                + tableNameMap.size());
        
        // Check the current version
        adapter = new MySQLAdaptor("localhost",
                FIConfiguration.getConfiguration().get("REACTOME_SOURCE_DB_NAME"),
                FIConfiguration.getConfiguration().get("DB_USER"),
                FIConfiguration.getConfiguration().get("DB_PWD")
                );
        conn = adapter.getConnection();
        stmt = conn.createStatement();
        rs = stmt.executeQuery("SHOW TABLES");
        while (rs.next()) {
            String tableName = rs.getString(1);
            String oldName = tableNameMap.get(tableName);
            if (oldName == null) {
                System.err.println("Cannot find table name for " + tableName);
            }
            else if (!oldName.equals(tableName)) {
                System.out.println(tableName + " -> " + oldName);  
                // Need a two-step approach 
                String tempName = oldName + "_temp_fix";
                String renameSql = String.format("ALTER TABLE `%s` RENAME TO `%s`;", tableName, tempName);
                System.out.println(renameSql);
                stmt.execute(renameSql);
                renameSql = String.format("ALTER TABLE `%s` RENAME TO `%s`;", tempName, oldName);
                System.out.println(renameSql);
                stmt.execute(renameSql);
            }
            else {
                System.out.println("Table name " + tableName + " is the same as the previous version.");
            }
        }
        rs.close();
        stmt.close();   
//        conn.close();
    }
    
    /**
     * Convert original INNODB FI database into myisam.
     * @throws Exception
     */
    @Test
    public void dumpFIDatabase() throws Exception {
        String fiDbName = FIConfiguration.getConfiguration().get("FI_DB_NAME");
        dump(fiDbName);
    }
    
    public void dumpReactomeSourceDatabase() throws Exception {
        String fiDbName = FIConfiguration.getConfiguration().get("REACTOME_SOURCE_DB_NAME");
        dump(fiDbName);
    }
    
    private void dump(String dbName) throws Exception {
        logger.info("Dump database " + dbName + "...");
        String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR") + File.separator + dbName + ".sql";
        String command = FIConfiguration.getConfiguration().get("MYSQLDUMP") + 
                         " -u" + FIConfiguration.getConfiguration().get("DB_USER") + 
                         " -p" + FIConfiguration.getConfiguration().get("DB_PWD") + 
                         " " + dbName;
        Process process = Runtime.getRuntime().exec(command);
        InputStream output = process.getInputStream();
        InputStreamReader isr = new InputStreamReader(output);
        BufferedReader bis = new BufferedReader(isr);
        String line = null;
        PrintWriter pr = new PrintWriter(fileName);
        while ((line = bis.readLine()) != null) {
            if (line.contains("ENGINE=InnoDB")) {
//                System.out.println(line);
                line = line.replaceAll("InnoDB", "MyISAM");
            }
            pr.println(line);
        }
        bis.close();
        isr.close();
        output.close();
        pr.close();
    }
    
    @Test
    public void dumpDatabases() throws Exception {
        dumpFIDatabase();
        dumpReactomeSourceDatabase();
    }
    
}
