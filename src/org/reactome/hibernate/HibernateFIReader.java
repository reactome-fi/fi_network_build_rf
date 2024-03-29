/*
 * Created on Oct 14, 2008
 *
 */
package org.reactome.hibernate;

import java.io.IOException;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.Statement;
import java.util.*;

import org.apache.log4j.Logger;
import org.hibernate.Query;
import org.hibernate.Session;
import org.hibernate.SessionFactory;
import org.junit.Test;
import org.reactome.data.ProteinIdFilters;
import org.reactome.data.UniProtAnalyzer;
import org.reactome.fi.FIFileAnalyzer;
import org.reactome.fi.ProteinAndInteractionCount;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.FileUtility;
import org.reactome.fi.util.InteractionUtilities;
import org.reactome.funcInt.Evidence;
import org.reactome.funcInt.Interaction;
import org.reactome.funcInt.Protein;
import org.reactome.funcInt.ReactomeSource;

/**
 * This class is used to query FIs based on the Hibernate APIs. This class is refactored from
 * HibernateFIAnalyzer.
 * @author wgm
 *
 */
public class HibernateFIReader extends HibernateFIPersistence {
    private static final Logger logger = Logger.getLogger(HibernateFIReader.class);
    
    private final Double CUT_OFF_VALUE = new Double(FIConfiguration.getConfiguration().get("CUT_OFF_VALUE"));

    public HibernateFIReader() {     
    }
    
    public void countInteractionsForProteins() throws Exception {
        Class.forName("com.mysql.jdbc.Driver");
        Connection connection = DriverManager.getConnection("jdbc:mysql://localhost:3306/FunctionalInteractions_v2",
                                                            "root",
                                                            "macmysql01");
        // Get the total of Proteins
        Statement stat = connection.createStatement();
        ResultSet results = stat.executeQuery("select max(dbId) from Protein");
        results.next();
        int maxProteinId = results.getInt(1);
        System.out.println("Maximum dbId of Protein: " + maxProteinId);
        results.close();
        stat.close();
        // These two for the total counting
        String query1 = "select count(*) from Interaction where firstProtein = ?";
        PreparedStatement prepStat1 = connection.prepareStatement(query1);
        String query2 = "select count(*) from Interaction where secondProtein = ?";
        PreparedStatement prepStat2 = connection.prepareStatement(query2);
        // These two are used for removing interactions with probabilities less than 0.6 (exclusive)
        String query3 = "select count(*) from Interaction i, Evidence e where i.firstProtein = ? " +
                        "and i.evidence = e.dbId and e.probability < 0.6";
        PreparedStatement prepStat3 = connection.prepareStatement(query3);
        String query4 = "select count(*) from Interaction i, Evidence e where i.secondProtein = ? " +
                        "and i.evidence = e.dbId and e.probability < 0.6";
        PreparedStatement prepStat4 = connection.prepareStatement(query4);
        PreparedStatement[] prepStat = new PreparedStatement[] {
                prepStat1,
                prepStat2,
                prepStat3,
                prepStat4
        };
        long time1 = System.currentTimeMillis();
        int counter = 0;
        // Output file
        String output = "results/ProteinInteractionCount0_6.txt";
        FileUtility fu = new FileUtility();
        fu.setOutput(output);
        //int c = 0;
        int value = 0;
        for (int i = 1; i < maxProteinId + 1; i++) {
            counter = 0;
            for (int j = 0; j < prepStat.length; j++) {
                // Total number
                prepStat[j].setLong(1, i);
                results = prepStat[j].executeQuery();
                results.next();
                value = results.getInt(1);
                results.close();
                if (j < 2)
                    counter += value;
                else
                    counter -= value;
            }
            fu.printLine(i + "\t" + counter);
            //if (c > 30)
            //    break;
            //c ++;
        }
        fu.close();
        long time2 = System.currentTimeMillis();
        System.out.println("Total Time: " + (time2 - time1));
        for (PreparedStatement stat1 : prepStat)
            stat1.close();
        connection.close();
    }
    
    /**
     * Query FIs based on gene names. Gene names are stored in Protein as shortNames.
     * @param name1
     * @param name2
     * @param session
     * @return
     * @throws Exception
     */
    public List<Interaction> queryFIsBasedOnGeneNames(String name1, 
                                                      String name2,
                                                      Session session) throws Exception {
        List<Interaction> rtn = new ArrayList<Interaction>();
        Query query = session.createQuery("FROM Interaction as i WHERE i.firstProtein.shortName = ? AND "
                                                                    + "i.secondProtein.shortName = ?");
        query.setString(0, name1);
        query.setString(1, name2);
        List interactions = query.list();
        if (interactions != null && interactions.size() > 0)
            rtn.addAll(interactions);
        // Do a reverse search
        query.setString(0, name2);
        query.setString(1, name1);
        interactions = query.list();
        if (interactions != null && interactions.size() > 0)
            rtn.addAll(interactions);
        return rtn;
    }
    
    @Test
    public void testQueryFIsBasedOnNamePair() throws Exception {
        String name1 = "TP53";
        String name2 = "BRCA1";
        String name3 = "BRCA2";
        String name4 = "RB1";
        initSession();
        Session session = sessionFactory.openSession();
        long time1 = System.currentTimeMillis();
        List<Interaction> rtn = new ArrayList<Interaction>();
        Query query = session.createQuery("SELECT p.dbId from Protein as p WHERE p.shortName in (?, ?)");
        query.setString(0, name1);
        query.setString(1, name2);
        List list = query.list();
        String listText = InteractionUtilities.joinStringElements(",", list);
        query = session.createQuery("FROM Interaction as i WHERE i.firstProtein in (" + listText + ") and i.secondProtein in (" + listText + ")");
        list = query.list();
        System.out.println("Total interactions: " + list.size());
        long time2 = System.currentTimeMillis();
        System.out.println("Dual query: " + (time2 - time1));
        
        query = session.createQuery("FROM Interaction as i WHERE i.firstProtein.shortName in ('TP53', 'BRCA1') and i.secondProtein.shortName in ('TP53', 'BRCA1')"); 
//        query.setString(0, name1);
//        query.setString(1, name2);
//        query.setString(2, name2);
//        query.setString(3, name1);
        List interactions = query.list();
        if (interactions != null && interactions.size() > 0)
            rtn.addAll(interactions);
//        // Do a reverse search
//        query.setString(0, name2);
//        query.setString(1, name1);
//        interactions = query.list();
//        if (interactions != null && interactions.size() > 0)
//            rtn.addAll(interactions);

        System.out.println("Total interctions: " + rtn.size());
        long time3 = System.currentTimeMillis();
        System.out.println("Time for query: " + (time3 - time2));
        session.close();
    }
    
    public void calculateRankFrequences() throws IOException {
        String fileName = "results/ProteinInteractionCount0_6.txt";
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = fu.readLine();
        int[] frequences = new int[800];
        int index = 0;
        String id, edges;
        int edgeNumber;
        while ((line = fu.readLine()) != null) {
            index = line.indexOf("\t");
            id = line.substring(0, index);
            edges = line.substring(index + 1);
            edgeNumber = Integer.parseInt(edges);
            frequences[edgeNumber] ++;
        }
        fu.close();
        // Get the total connections
        int total = 0;
        for (int tmp : frequences)
            total += tmp;
        System.out.println("Total: " + total);
        String output = "results/ConnectionFrequences0_6.txt";
        fu.setOutput(output);
        fu.printLine("edgeNumber\tfrequence");
        // Don't want to include 0: meaningless
        for (int i = 1; i < frequences.length; i++) {
            if (frequences[i] == 0)
                continue;
            fu.printLine(i + "\t" + (double)frequences[i] / total);
        }
        fu.close();
        // Want to sort the frequences
        // to exclude the sort for the first element
        frequences[0] = 0;
        Arrays.sort(frequences);
        output = "results/RankOfFrequences0_6.txt";
        fu.setOutput(output);
        fu.printLine("rankOfFrequence\tfrequence");
        int rank = 1;
        for (int i = frequences.length - 1; i > -1; i--) {
            if (frequences[i] == 0)
                break; // the last one
            fu.printLine(rank + "\t" + (double)frequences[i] / total);
            rank ++;
        }
        fu.close();
    }
    
    /**
     * This method is used to generate FIs from a hibernated FI database.
     * @throws Exception
     */
    public void generatePathwayFIFileInHibernate() throws Exception {
        initSession();
        Session session = sessionFactory.openSession();
        Query query = session.createQuery("FROM Interaction as i WHERE i.evidence is null");
        List list = query.list();
        System.out.println("Total interactions from pathways: " + list.size());
        Set<String> fis = new HashSet<String>();
        for (Iterator it = list.iterator(); it.hasNext();) {
            Interaction interaction = (Interaction) it.next();
            Protein protein1 = interaction.getFirstProtein();
            Protein protein2 = interaction.getSecondProtein();
            fis.add(protein1.getPrimaryAccession() + " " + protein2.getPrimaryAccession());
        }
        session.close();
        String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR") + "FI73_Pathway_From_DB_041808.txt";
        new FIFileAnalyzer().saveFIInOrder(fis, fileName);
    }
    
    protected String getFIInName(Interaction interaction) {
        Protein protein1 = interaction.getFirstProtein();
        Protein protein2 = interaction.getSecondProtein();
        //String name1 = protein1.getLabel();
        //String name2 = protein2.getLabel();
        // As of April 17, 2009, use geneName only
        String name1 = protein1.getShortName();
        String name2 = protein2.getShortName();
        if (name1 == null || name2 == null) {
            // There are too many cases like this. Just don't export them!
//            logger.warn("Interaction " + interaction.getDbId()  + " cannot be converted to FI in Gene: " + 
//                         protein1.getPrimaryAccession() + " - " + protein2.getPrimaryAccession());
            // Do some check
            return null;
        }
        // Note: Make sure all names are in upper case!!!
        name1 = name1.toUpperCase();
        name2 = name2.toUpperCase();
        int compare = name1.compareTo(name2);
        String fi = null;
        // Need to use tab for delimited in case space is used
        // in names
        if (compare < 0)
            fi = name1 + "\t" + name2;
        else if (compare > 0)
            fi = name2 + "\t" + name1;
        return fi;
    }
    
    /**
     * Dump a map from FIs in genes to their sources in DB_IDs.
     * @throws Exception
     */
    @Test
    public void generateFIInGeneSourceFile() throws Exception {
        initSession();
        Session session = sessionFactory.openSession();
        long time1 = System.currentTimeMillis();
        Query query = session.createQuery("FROM Interaction as i WHERE i.evidence is null");
        List<?> list = query.list();
        logger.info("Total interactions from pathways: " + list.size());
        Interaction interaction = null;
        // Use String to use other methods
        Map<String, Set<String>> fiToSources = new HashMap<String, Set<String>>();
        for (Iterator<?> it = list.iterator(); it.hasNext();) {
            interaction = (Interaction) it.next();
            String fi = getFIInName(interaction);
            if (fi == null)
                continue;
            Set<ReactomeSource> fiSources = interaction.getReactomeSources();
            if (fiSources == null)
                continue;
            for (ReactomeSource fiSource : fiSources) {
                //                if (fiSource.getDataSource().equals("Reactome"))
                InteractionUtilities.addElementToSet(fiToSources, fi, fiSource.getReactomeId() + "");
            }
        }
        long time2 = System.currentTimeMillis();
        logger.info("Time for getting interactions: " + (time2 - time1));
        session.close();
        logger.info("Total Pathway FIs in Genes with sources: " + fiToSources.size());
//        logger.info("Total Reactome Pathway FIs in Genes with sources: " + fiToSources.size());
        FileUtility fu = new FileUtility();
        fu.saveSetMap(fiToSources, 
                      FIConfiguration.getConfiguration().get("GENE_FI_PATHWAY_SOURCES_FILE_NAME"));
    }
    
    protected List fetchPredictedFIs(Session session) throws Exception {
        Query query = session.createQuery("FROM Interaction as i WHERE i.evidence.score >= ?");
        query.setDouble(0, CUT_OFF_VALUE);
        List list = query.list();
        System.out.println("Total interactions from prediction: " + list.size());
        return list;
    }
    
    protected List fetchAnnotatedFIs(Session session) throws Exception {
        Query query = session.createQuery("FROM Interaction as i WHERE i.evidence is null");
        List list = query.list();
        System.out.println("Total interactions from pathways: " + list.size());
        return list;
    }
    
    /**
     * This method is used to generate three interaction files in gene names using
     * hibernate APIs.
     * @throws Exception
     */
    @Test
    public void generateFIFileInGeneInHibernate() throws Exception {
        initSession();
        Session session = sessionFactory.openSession();
        long time1 = System.currentTimeMillis();
        List list = fetchPredictedFIs(session);
        Interaction interaction = null;
        //String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR + "FI73_041408.txt";
        // Used to hold FIs to count how many proteins and interactions
        Set<String> predictedFIs = new HashSet<String>();
        Set<String> pathwayFIs = new HashSet<String>();
        Set<String> allFIs = new HashSet<String>();
        // Start counting and output
        for (Iterator it = list.iterator(); it.hasNext();) {
            interaction = (Interaction) it.next();
            String fi = getFIInName(interaction);
            if (fi == null)
                continue;
            predictedFIs.add(fi);
            allFIs.add(fi);
        }
        list = fetchAnnotatedFIs(session);
        for (Iterator it = list.iterator(); it.hasNext();) {
            interaction = (Interaction) it.next();
//            // Avoid FIs from panther and INOH
//            Set<ReactomeSource> sources = interaction.getReactomeSources();
//            boolean shouldScape = true;
//            for (ReactomeSource src : sources) {
//                if (!src.getDataSource().equals("pantherdb") &&
//                    !src.getDataSource().equals("INOH")) {
//                    shouldScape = false;
//                    break;
//                }
//            }
//            if (shouldScape)
//                continue;
            String fi = getFIInName(interaction);
            if (fi != null) {
                allFIs.add(fi);
                pathwayFIs.add(fi);
            }
        }
        long time2 = System.currentTimeMillis();
        System.out.println("Time for getting interactions: " + (time2 - time1));
        session.close();
//        System.out.println("Before name normalization...");
        System.out.println("Total predicted FIs: " + predictedFIs.size());
        System.out.println("Total pathway FIs:" + pathwayFIs.size());
        System.out.println("Total FIs: " + allFIs.size());
//        EntrezGeneAnalyzer entrezAnalyzer = new EntrezGeneAnalyzer();
//        allFIs = entrezAnalyzer.normalizeFIInNames(allFIs);
//        predicatedFIs = entrezAnalyzer.normalizeFIInNames(predicatedFIs);
//        pathwayFIs = entrezAnalyzer.normalizeFIInNames(pathwayFIs);
//        // Calculate numbers
//        // FIs
//        System.out.println("After name normalization...");
//        System.out.println("Total predicated FIs: " + predicatedFIs.size());
//        System.out.println("Total pathway FIs:" + pathwayFIs.size());
//        System.out.println("Total FIs: " + allFIs.size());
        // Proteins in FIs
        Set<String> predicatedProteins = InteractionUtilities.grepIDsFromInteractions(predictedFIs);
        Set<String> pathwayProteins = InteractionUtilities.grepIDsFromInteractions(pathwayFIs);
        Set<String> allProteins = InteractionUtilities.grepIDsFromInteractions(allFIs);
        System.out.println("Total predicted proteins: " + predicatedProteins.size());
        System.out.println("Total pathway proteins: " + pathwayProteins.size());
        System.out.println("Total proteins: " + allProteins.size());
        // Output FI files
        FIFileAnalyzer fiFileAnalyzer = new FIFileAnalyzer();
        fiFileAnalyzer.saveFIInOrder(allFIs, 
                                     FIConfiguration.getConfiguration().get("GENE_FI_FILE_NAME"));
        fiFileAnalyzer.saveFIInOrder(predictedFIs,
                                     FIConfiguration.getConfiguration().get("GENE_FI_PREDICTED_FILE_NAME"));
        fiFileAnalyzer.saveFIInOrder(pathwayFIs,
                                     FIConfiguration.getConfiguration().get("GENE_FI_PATHWAY_FILE_NAME"));
//         All files
//        fiFileAnalyzer.saveFIInOrder(allFIs, 
//                                     FIConfiguration.getConfiguration().get("RESULT_DIR + "FI73InGene_102908.txt");
//        fiFileAnalyzer.saveFIInOrder(allFIs, 
//                                     FIConfiguration.getConfiguration().get("RESULT_DIR + "FI73InGene_062008_NO_P_I.txt");
//        fiFileAnalyzer.saveFIInOrder(predictedFIs,
//                                     FIConfiguration.getConfiguration().get("RESULT_DIR + "FI73InGene_Predicated_062008_NO_P_I.txt");
//        fiFileAnalyzer.saveFIInOrder(pathwayFIs, 
//                                     FIConfiguration.getConfiguration().get("RESULT_DIR + "FI73InGene_Pathway_062008_NO_P_I.txt");
    }
    
    /**
     * Query for Interactions for a passed gene or protein name.
     * @param name
     * @param session
     * @return
     * @throws Exception
     */
    public List<Interaction> queryFIsForName(String name, 
                                             Session session) throws Exception {
        List<Interaction> rtn = new ArrayList<Interaction>();
        // Get protein first
        Query query = session.createQuery("FROM Interaction as i WHERE i.firstProtein.shortName = ?");
        query.setString(0, name);
        List list = query.list();
        if (list != null) {
            for (Iterator it = list.iterator(); it.hasNext();) {
                rtn.add((Interaction)it.next());
            }
        }
        query = session.createQuery("FROM Interaction as i WHERE i.secondProtein.shortName = ?");
        query.setString(0, name);
        list = query.list();
        if (list != null) {
            for (Iterator it = list.iterator(); it.hasNext();) {
                rtn.add((Interaction)it.next());
            }
        }
        return rtn;
    }
    
    @Test
    public void testQueryFIsForAccession() throws Exception {
        String accession = "Q8N726";
        initSession();
        Session session = sessionFactory.openSession();
        List<Interaction> fis = queryFIsForAccessio(accession, session);
        System.out.println("Total FIs for " + accession + ": " + fis.size());
        for (Interaction fi : fis) {
            Evidence evidence = fi.getEvidence();
//            System.out.println(fi.getFirstProtein().getPrimaryAccession() + "\t" + 
//                               fi.getSecondProtein().getPrimaryAccession() + "\t" + 
//                               fi.getFirstProtein().getShortName() + "\t" + 
//                               fi.getSecondProtein().getShortName() + "\t" + 
//                               (evidence == null ? "extractedFI" : 
//                                                   ("predictedFI" + "\t" + evidence.getProbability() + "\t" + 
//                                                                           evidence.getHumanInteraction() + "\t" + 
//                                                                           evidence.getMousePPI() + "\t" + 
//                                                                           evidence.getDmePPI() + "\t" + 
//                                                                           evidence.getScePPI() + "\t" + 
//                                                                           evidence.getCelPPI() + "\t" +
//                                                                           evidence.getPfamDomainInt())));
            if (evidence != null)
                continue;
            Set<ReactomeSource> sources = fi.getReactomeSources();
            for (ReactomeSource source : sources) {
                if (source.getDataSource().equals("Reactome")) {
                    System.out.println(fi.getFirstProtein().getPrimaryAccession() + "\t" + 
                                       fi.getSecondProtein().getPrimaryAccession() + "\t" + 
                                       fi.getFirstProtein().getShortName() + "\t" + 
                                       fi.getSecondProtein().getShortName());
                    break;
                }
            }
        }
        session.close();
    }
    
    public List<Interaction> queryFIsForAccessio(String accession, 
                                                 Session session) throws Exception {
        List<Interaction> rtn = new ArrayList<Interaction>();
        // Get protein first
        Query query = session.createQuery("FROM Interaction as i WHERE i.firstProtein.primaryDbReference.accession = ?");
        query.setString(0, accession);
        List list = query.list();
        if (list != null) {
            for (Iterator it = list.iterator(); it.hasNext();) {
                rtn.add((Interaction)it.next());
            }
        }
        query = session.createQuery("FROM Interaction as i WHERE i.secondProtein.primaryDbReference.accession = ?");
        query.setString(0, accession);
        list = query.list();
        if (list != null) {
            for (Iterator it = list.iterator(); it.hasNext();) {
                rtn.add((Interaction)it.next());
            }
        }
        return rtn;
    }
    
    /**
     * This method is used to generate protein names to UniProt accession numbers.
     * @throws Exception
     */
    public Map<String, Set<String>> generateProteinNameToAccession() throws Exception {
        initSession();
        Session session = sessionFactory.openSession();
        Query query = session.createQuery("FROM Protein");
        List list = query.list();
        Map<String, Set<String>> nameToIds = new HashMap<String, Set<String>>();
        for (Iterator it = list.iterator(); it.hasNext();) {
            Protein protein = (Protein) it.next();
            String accession = protein.getPrimaryAccession();
            String name = protein.getLabel();
            Set<String> ids = nameToIds.get(name);
            if (ids == null) {
                ids = new HashSet<String>();
                nameToIds.put(name, ids);
            }
            ids.add(accession);
        }
        session.close();
        return nameToIds;
    }
    
    public Map<String, String> generateAccessionToProteinNames() throws Exception {
        initSession();
        Session session = sessionFactory.openSession();
        Query query = session.createQuery("FROM Protein");
        List list = query.list();
        Map<String, String> idToName = new HashMap<String, String>();
        for (Iterator it = list.iterator(); it.hasNext();) {
            Protein protein = (Protein) it.next();
            String id = protein.getPrimaryAccession();
            String name = protein.getLabel();
            idToName.put(id, name);
        }
        session.close();
        return idToName;
    }
    
    @Test
    public void generateAccessionToProteinNameMap() throws Exception {
        Map<String, String> idToName = generateAccessionToProteinNames();
//        String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR") + "ProteinAccessionToName_070110.txt";
        String fileName = FIConfiguration.getConfiguration().get("PROTEIN_ACCESSION_TO_NAME_FILE");
        new FileUtility().exportMap(idToName, fileName);
    }
    
    /**
     * This method is used to analyze the mapping from name to ids.
     * @throws Exception
     */
    @Test
    public void analyzeProteinNameToAccessions() throws Exception {
        Map<String, Set<String>> idToNames = generateProteinNameToAccession();
        // Want to reverse the mapping
        final Map<String, Set<String>> nameToIds = InteractionUtilities.switchKeyValues(idToNames);
        // Do a sorting
        List<String> names = new ArrayList<String>(nameToIds.keySet());
        Collections.sort(names, new Comparator<String>() {
            public int compare(String name1,
                               String name2) {
                Set<String> ids1 = nameToIds.get(name1);
                Set<String> ids2 = nameToIds.get(name2);
                return ids2.size() - ids1.size(); 
            }
        });
        for (String name : names) {
            Set<String> ids = nameToIds.get(name);
            if (ids.size() == 1)
                continue;
            System.out.println(name + ": " + ids.size());
        }
    }
    
    /**
     * Generate an FI file for the FI database prior version 4
     * @throws Exception
     */
    public void generateInteractionFile() throws Exception {
        // Used to check UniProt Ids
        UniProtAnalyzer uniAnalyzer = new UniProtAnalyzer();
        Map<String, String> uniProtIds = uniAnalyzer.loadUniProtIDsMap();
        Set<String> uniSet = uniProtIds.keySet();
        
        Class.forName("com.mysql.jdbc.Driver");
        Connection connection = DriverManager.getConnection("jdbc:mysql://localhost:3306/FunctionalInteractions_v2",
                                                            "root",
                                                            "macmysql01");
        String query = "SELECT p1.accession, p2.accession From Interaction i, Protein p1, Protein p2 " +
                "where i.firstProtein = p1.dbId and i.secondProtein = p2.dbId AND p1.accession != -1 " +
                "AND i.evidence is null order by i.dbId";
        Set<String> pathwayInteractions = new HashSet<String>();
        FileUtility fu = new FileUtility();
        //fu.setOutput("results/v2/FIInteractions06.sif");
        // Before update: regulators are not included in pathway participants
        // fu.setOutput("results/v2/FIInteractions07.txt");
        // Regulators are includes
        //fu.setOutput("results/v2/FIInteractions67.txt");
        //fu.setOutput("results/v2/FIInteractions57.txt");
        // the cutoff values should be 0.60 as of Dec 17, 2007
        double cutoff = 0.73;
        Statement stat = connection.createStatement();
        ResultSet resultSet = stat.executeQuery(query);
        while (resultSet.next()) {
            String acc1 = resultSet.getString(1);
            if (!uniAnalyzer.isHumanID(uniSet, acc1))
                continue;
            String acc2 = resultSet.getString(2);
            if (!uniAnalyzer.isHumanID(uniSet, acc2))
                continue;
            // Need to make sure the first ids are used
            if (uniProtIds.containsKey(acc1))
                acc1 = uniProtIds.get(acc1);
            if (uniProtIds.containsKey(acc2))
                acc2 = uniProtIds.get(acc2);
            int compare = acc1.compareTo(acc2);
            if (compare < 0)
                pathwayInteractions.add(acc1 + " " + acc2);
            else if (compare > 0)
                pathwayInteractions.add(acc2 + " " + acc1);
        }
        resultSet.close();
        query = "SELECT p1.accession, p2.accession From Interaction i, Evidence e, Protein p1, Protein p2 " +
        "where i.firstProtein = p1.dbId and i.secondProtein = p2.dbId AND p1.accession != -1 " +
        "AND i.evidence = e.dbId AND e.probability >= " + cutoff + " order by i.dbId";
        Set<String> ppiInteractions = new HashSet<String>();
        resultSet = stat.executeQuery(query);
        while (resultSet.next()) {
            String acc1 = resultSet.getString(1);
            if (!uniAnalyzer.isHumanID(uniSet, acc1))
                continue;
            String acc2 = resultSet.getString(2);
            if (!uniAnalyzer.isHumanID(uniSet, acc2))
                continue;
            // Need to make sure the first ids are used
            if (uniProtIds.containsKey(acc1))
                acc1 = uniProtIds.get(acc1);
            if (uniProtIds.containsKey(acc2))
                acc2 = uniProtIds.get(acc2);
            int compare = acc1.compareTo(acc2);
            if (compare < 0)
                ppiInteractions.add(acc1 + " " + acc2);
            else if (compare > 0)
                ppiInteractions.add(acc2 + " " + acc1);
        }
        resultSet.close();
        stat.close();
        connection.close();
        // This is used to output FIs from pathways only
        String pathwayFileName = FIConfiguration.getConfiguration().get("RESULT_DIR") + "FIInteractions73_021108_Pathway.txt";
        FileUtility fu1 = new FileUtility();
        fu1.setOutput(pathwayFileName);
        List<String> pathwaysInList = new ArrayList<String>(pathwayInteractions);
        Collections.sort(pathwaysInList);
        for (String line : pathwaysInList)
            fu1.printLine(line);
        fu1.close();
        List<String> list = new ArrayList<String>(ppiInteractions);
        Collections.sort(list);
        // This is used to ouput FIs from Human PPIs
        String ppiFileName = FIConfiguration.getConfiguration().get("RESULT_DIR") + "FIInteractions73_021108_PPI.txt";
        fu1.setOutput(ppiFileName);
        for (String line : list)
            fu1.printLine(line);
        fu1.close();
        // Output all
        Set<String> interactions = new HashSet<String>();
        interactions.addAll(ppiInteractions);
        interactions.addAll(pathwayInteractions);
        list = new ArrayList<String>(interactions);
        String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR") + "FIInteractions73_021108.txt";
        fu1.setOutput(fileName);
        for (String line : list)
            fu1.printLine(line);
        fu1.close();
        System.out.println("Size: pathways: " + pathwayInteractions.size() +
                           ", ppi: " + ppiInteractions.size() + ", " +
                           ", total: " + interactions.size());
                           
    }
    
    private List<String> fetchAccessionsForProteinIds(Set<Long> ids,
                                                      Connection connection,
                                                      boolean isOldDb) throws Exception {
        String query = null;
        if (isOldDb) {
            query = "SELECT accession from protein where dbId IN ";
        }
        else {
            query = "select accession from protein p, dbreference d where p.primaryDbReference = d.dbId and p.dbId IN ";
        }
        // Need to consider a long ids
        StringBuilder builder = new StringBuilder();
        builder.append("(");
        for (Iterator<Long> it = ids.iterator(); it.hasNext();) {
            Long id = it.next();
            builder.append(id);
            if (it.hasNext())
                builder.append(",");
        }
        builder.append(")");
        query = query + builder.toString();
        Statement stat = connection.createStatement();
        ResultSet resultset = stat.executeQuery(query);
        List<String> accessions = new ArrayList<String>();
        while (resultset.next())
            accessions.add(resultset.getString(1));
        resultset.close();
        stat.close();
        return accessions;
    }
    
    /**
     * This method is used to count proteins and interactions from a FI database
     * by using JDBC directly.
     */
    public void countProteinsAndInteractionsDirectly() throws Exception {
        Class.forName("com.mysql.jdbc.Driver");
        String url = null;
        String[] dbNames = new String[] {
                "FunctionalInteractions_v3",
                "FunctionalInteractions_v4",
                "FI_v4_Pathway_PPI"
        };
        UniProtAnalyzer uniAnalyzer = new UniProtAnalyzer();
        Set<String> swissProtIds = uniAnalyzer.loadSwissProtIds();
        for (String dbName : dbNames) {
            // FunctionalInteraction database version 3
            url = "jdbc:mysql://localhost:3306/" + dbName;
            // FunctionalInteraction database version 4: pathway FIs and predicated FIs from human PPI
            //url = "jdbc:mysql://localhost:3306/FI_v4_Pathway_PPI";
            // FunctionalInteraction database version 4: FIs from pathway, human PPIs and other sources
            //url = "jdbc:mysql://localhost:3306/FunctionalInteractions_v4";
            Connection connection = DriverManager.getConnection(url,
                                                                "root",
                                                                "macmysql01");
            Set<Long> proteins = new HashSet<Long>();
            int interactionCount = 0;
            // This query is not extracted pathway FIs
            String query = "SELECT firstProtein, secondProtein FROM interaction where evidence is null";
            Statement stat = connection.createStatement();
            ResultSet resultSet = stat.executeQuery(query);
            while (resultSet.next()) {
                Long protein1 = resultSet.getLong(1);
                Long protein2 = resultSet.getLong(2);
                proteins.add(protein1);
                proteins.add(protein2);
                interactionCount ++;
            }
            resultSet.close();
            System.out.println("Database: " + dbName);
            System.out.println("Total proteins from pathways: " + proteins.size());
            System.out.println("Total interactions from pathways: " + interactionCount);
            boolean isOld = dbName.endsWith("_v3");
            List<String> accessions = fetchAccessionsForProteinIds(proteins, connection, isOld);
            accessions.retainAll(swissProtIds);
            System.out.println("Total SwissProt ids: " + accessions.size() + " (" + accessions.size() / (double) swissProtIds.size() + ")");
            // This query is for predicated FIs
            query = "SELECT i.firstProtein, i.secondProtein FROM interaction i, evidence e where i.evidence = e.dbId and e.probability >= 0.73";
            resultSet = stat.executeQuery(query);
            while (resultSet.next()) {
                Long protein1 = resultSet.getLong(1);
                Long protein2 = resultSet.getLong(2);
                proteins.add(protein1);
                proteins.add(protein2);
                interactionCount ++;
            }
            resultSet.close();
            stat.close();
            System.out.println("Total proteins: " + proteins.size());
            System.out.println("Total interactions: " + interactionCount);
            accessions = fetchAccessionsForProteinIds(proteins, connection, isOld);
            accessions.retainAll(swissProtIds);
            System.out.println("Total SwissProt ids: " + accessions.size() + " (" + accessions.size() / (double) swissProtIds.size() + ")");
            System.out.println();
        }
    }
    
    /**
     * This method should be in the top-level for testing. However, for some reason it cannot work.
     * It must be in the subclass!
     * @throws Exception
     */
    @Test
    public void testSetting() throws Exception {
        System.out.println("Starting checking...");
        initSession();
        Session session = sessionFactory.openSession();
        session.close();
    }
    
    @Test
    public void checkFIsForGene() throws Exception {
        initSession();
        Session session = sessionFactory.openSession();
        //String name = "TARDBP";
        String name = "AHNAK";
        List<Interaction> interactions = queryFIsForName(name, session);
        System.out.println("FIs for " + name + ": " + interactions.size());
        for (Interaction i : interactions) {
            Protein p1 = i.getFirstProtein();
            Protein p2 = i.getSecondProtein();
            System.out.println(p1.getLabel() + "\t" + p2.getLabel() + "\t" + i.getEvidence());
        }
        session.close();
    }
    
//    @Test
//    public void checkProteins() throws Exception {
//        initSession();
//        Session session = sessionFactory.openSession();
//        Query query = session.createQuery("FROM Protein");
//        List list = query.list();
//        Set<String> dbAccs = new HashSet<String>();
//        for (Iterator it = list.iterator(); it.hasNext();) {
//            Protein protein = (Protein) it.next();
//            dbAccs.add(protein.getPrimaryAccession());
//        }
//        session.close();
//        System.out.println("Total db accessions: " + dbAccs.size());
//        // Get proteins from the local files
//        Set<String> localFIs = new FIFileAnalyzer().loadFIs();
//        Set<String> localAccs = InteractionUtilities.grepIDsFromInteractions(localFIs);
//        System.out.println("Total local accessions: " + localAccs.size());
//        Set<String> localCopy = new HashSet<String>(localAccs);
//        localCopy.removeAll(dbAccs);
//        System.out.println("Accession in local but not in db: " + localCopy);
//        dbAccs.removeAll(localAccs);
//        System.out.println("Accession in db but not in local: " + dbAccs);
//    }
    
    @Test
    public void checkEvidences() throws Exception {
        initSession();
        Session session = sessionFactory.openSession();
        Query query = session.createQuery("FROM Evidence");
        List list = query.list();
        Set<String> dbAccs = new HashSet<String>();
        for (Iterator it = list.iterator(); it.hasNext();) {
            Evidence evidence = (Evidence) it.next();
            System.out.println(evidence.getDbId() + ", humanInteraction: " + evidence.getHumanPPI() + ", " + evidence.getScore());
        }
        session.close();
    }
    
    @Test
    public void checkPantherInteractions() throws Exception {
        initSession();
        Session session = sessionFactory.openSession();
        // Check how many unique FIs extracted from Panther
        Query query = session.createQuery("FROM Interaction");
        @SuppressWarnings("unchecked")
        List<Interaction> interactions = query.list();
        int total = 0;
        for (Interaction interaction : interactions) {
            Set<ReactomeSource> sources = interaction.getReactomeSources();
            if (sources == null || sources.size() == 0)
                continue;
            if (sources.size() == 1) {
                ReactomeSource source = sources.iterator().next();
                if (source.getDataSource().equals("pantherdb"))
                    total ++;
            }
        }
        session.close();
        System.out.println("Interactions extracted only from Panther: " + total);
    }
    
    @Test
    public void checkTREDFIs() throws Exception {
        initSession();
        Session session = sessionFactory.openSession();
        Query query = session.createQuery("FROM Interaction");
        List<?> interactions = query.list();
        List<Interaction> tredInteractions = new ArrayList<Interaction>();
        for (Object obj : interactions) {
            Interaction interaction = (Interaction) obj;
            Set<ReactomeSource> sources = interaction.getReactomeSources();
            for (ReactomeSource source : sources){
                if (source.getDataSource().equals("TRED")) {
                    tredInteractions.add(interaction);
                    break;
                }
            }
        }
        System.out.println("Total interactions from TRED: " + tredInteractions.size());
        Set<String> fis = new HashSet<String>();
        for (Interaction interaction : tredInteractions) {
            Protein protein1 = interaction.getFirstProtein();
            Protein protein2 = interaction.getSecondProtein();
            System.out.println(protein1.getPrimaryAccession() + "\t" + protein2.getPrimaryAccession());
            fis.add(protein1.getPrimaryAccession() + "\t" + protein2.getPrimaryAccession());
        }
        ProteinIdFilters filters = new ProteinIdFilters();
        Set<String> normalized = filters.normalizeProteinPairs(fis);
        System.out.println("After filtering: " + normalized.size());
        session.close();
    }
    
    /**
     * Load gene name to protein name mapping.
     * @param geneNames
     * @return
     * @throws Exception
     */
    public Map<String, String> loadShortNameToNameMap(Collection<String> geneNames) throws Exception {
        Map<String, String> map = new HashMap<String, String>();
        initSession();
        Session session = sessionFactory.openSession();
        String queryString = "FROM " + Protein.class.getName() + " p WHERE p.shortName = ?";
        Query query = session.createQuery(queryString);
        for (String geneName : geneNames) {
            query.setParameter(0, geneName);
            List list = query.list();
            if (list == null || list.size() == 0)
                continue;
            Protein protein = (Protein) list.get(0);
            map.put(geneName, protein.getName());
        }
        session.close();
        return map;
    }
}
