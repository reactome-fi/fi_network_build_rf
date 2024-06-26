/*
 * Created on Oct 14, 2008
 *
 */
package org.reactome.fi;

import java.io.File;
import java.lang.reflect.Method;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;
import org.hibernate.Query;
import org.hibernate.Session;
import org.hibernate.Transaction;
import org.junit.Test;
import org.reactome.data.EncodeInteractionAnalyzer;
import org.reactome.data.ReactomeAnalyzer;
import org.reactome.data.ReactomeFuncInteractionExtractor;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.HibernateUtil;
import org.reactome.fi.util.PositiveChecker;
import org.reactome.fi.util.Value;
import org.reactome.funcInt.Evidence;
import org.reactome.funcInt.Interaction;
import org.reactome.funcInt.Protein;
import org.reactome.funcInt.ReactomeSource;
import org.reactome.hibernate.HibernateFIPersistence;
import org.reactome.weka.NaiveBayesClassifier;
import org.reactome.weka.WEKADataAnalyzer;
import org.reactome.weka.WEKAResultAnalyzer;

import weka.core.Instances;

/**
 * This class is used to write FIs to the database based on the Hibernate API. This class 
 * was refactored from HibernateFIAnalyzer.
 * @author wgm
 *
 */
public class FIDBBuilder extends HibernateFIPersistence {
    private static final Logger logger = Logger.getLogger(FIDBBuilder.class);
    private final double CUT_OFF_VALUE = Double.parseDouble(FIConfiguration.getConfiguration().get("CUT_OFF_VALUE"));
    
    public FIDBBuilder() {   
    }
    
    /**
     * Use this method to create the database schema for the FI network.
     * @throws Exception
     */
    @Test
    public void generateSchema() throws Exception {
        String configFileName = "resources/funcIntHibernate.cfg.xml";
        File configFile = new File(configFileName);
        sessionFactory = HibernateUtil.getSessionFactoryWithCreate(configFile);
        Session session = sessionFactory.openSession();
        Query query = session.createQuery("FROM " + Interaction.class.getName());
        List<?> list = query.list();
        logger.info("Interaction should be empty: " + list.size());
        if (list.size() > 0)
            throw new IllegalStateException("FI Database is not empty!");
        session.close();
    }
    
    @Test
    public void testDump() throws Exception {
        ReactomeFuncInteractionExtractor extractor = new ReactomeFuncInteractionExtractor();
        List<ReactomeAnalyzer> analyzerList = ReactomeAnalyzer.getPathwayDBAnalyzers();
        for (ReactomeAnalyzer a : analyzerList) {
            // Just want to check encode
            if (!(a instanceof EncodeInteractionAnalyzer))
                continue;
            extractor.setReactomeAnalyzer(a);
            extractor.extractFuncInteractions();
        }
        List<Interaction> interactions = extractor.getExtractedInteractions();
        logger.info("total interactions: " + interactions.size());
    }
    
    /**
     * Method to dump all interactions extracted from the databases into the FI database.
     * @throws Exception
     */
    @Test
    public void dump() throws Exception {
        long time1 = System.currentTimeMillis();
        ReactomeFuncInteractionExtractor extractor = new ReactomeFuncInteractionExtractor();
        List<ReactomeAnalyzer> analyzerList = ReactomeAnalyzer.getPathwayDBAnalyzers();
        for (ReactomeAnalyzer a : analyzerList) {
            if (a.getDataSource() == null)
                logger.info("Extract interactions from Reactome...");
            else
                logger.info("Extract interactions from " + a.getDataSource().getDisplayName() + "...");
            extractor.setReactomeAnalyzer(a);
            extractor.extractFuncInteractions();
        }
        List<Interaction> interactions = extractor.getExtractedInteractions();
        logger.info("total interactions: " + interactions.size());
        long time2 = System.currentTimeMillis();
        logger.info("Total time to extract: " + (time2 - time1));
        initSession();
        Session session = sessionFactory.getCurrentSession();
        Transaction tx = null;
        try {
            tx = session.beginTransaction();
            for (Interaction interaction : interactions) {
                session.persist(interaction);
            }
            tx.commit();
            if (session.isOpen())
                session.close();
        }
        catch(Exception e) {
            if (tx != null)
                tx.rollback();
            if (session.isOpen())
                session.close();
            throw e; // rethrow the exception to be caught by the caller.
        }
        long time3 = System.currentTimeMillis();
        logger.info("Total time for saving: " + (time3 - time2));
    }
    
    /**
     * This method should be invoked after the dump method.
     * @throws Exception
     */
    public void addEvidences() throws Exception {
        initSession();
        Session session = sessionFactory.openSession();
        Transaction tx = session.beginTransaction();
        Query query = session.createQuery("FROM " + Interaction.class.getName());
        int start = 0;
        int max = 5000;
        // Use pagination query to minimize the memory usage.
        query.setFirstResult(start);
        query.setMaxResults(max);
        long time1 = System.currentTimeMillis();
        List list = query.list();
        Interaction interaction = null;
        List<Interaction> intNeedEvidences = new ArrayList<Interaction>();
        boolean needEvidence = true;
        while (list.size() > 0) {
            //intNeedEvidences.clear();
            for (Iterator it = list.iterator(); it.hasNext();) {
                interaction = (Interaction) it.next();
                needEvidence = true;
                Set<ReactomeSource> sources = interaction.getReactomeSources();
                for (ReactomeSource source : sources) {
                    String dbName = source.getDataSource();
                    if (!dbName.equals("BIND") &&
                        !dbName.equals("IntAct") &&
                        !dbName.equals("HPRD")) {
                        needEvidence = false;
                        break;
                    }
                }
                if (needEvidence)
                    intNeedEvidences.add(interaction);
            }
            //session.clear(); // To empty cache to keep the memory usage small.
            start += list.size();
            query.setFirstResult(start);
            list = query.list();
            logger.info("Current interactions: " + start);
            //break;
        }
        long time2 = System.currentTimeMillis();
        logger.info("Time for getting interactions: " + (time2 - time1));
        logger.info("Instances needing evidence: " + intNeedEvidences.size());
        // Need to make change here just before session.clear. Otherwise,
        // changed states cannot save fast enough.
        addEvidences(intNeedEvidences);
        //addEvidenceForManuel(intNeedEvidences);
        long time3 = System.currentTimeMillis();
        logger.info("Time for adding evidences: " + (time3 - time2));
        int c = 0;
        for (Interaction tmp : intNeedEvidences) {
            session.persist(tmp.getEvidence());
//            c ++;
//            if (c % 500 == 0) {
//                session.flush();
//                session.clear();
//            }
        }
        session.flush();
        session.clear();
        long time4 = System.currentTimeMillis();
        logger.info("Time for saving evidences: " + (time4 - time3));
        // Update Interactions
        updateInteractionsForEvidences(intNeedEvidences, session);
        long time5 = System.currentTimeMillis();
        logger.info("Time for updating interactions: " + (time5 - time4));
        tx.commit();
        session.close();
    }
    
    private void updateInteractionsForEvidences(List<Interaction> interactions,
                                                Session session) throws Exception {
        Connection connection = session.connection();
        PreparedStatement stat = connection.prepareStatement("UPDATE Interaction SET evidence = ? WHERE dbId = ?");
        for (Interaction inter : interactions) {
            stat.setLong(1, inter.getEvidence().getDbId());
            stat.setLong(2, inter.getDbId());
            stat.executeUpdate();
        }
    }
    
    /**
     * This helper method is used to load all proteins in the database into a map, which
     * is keyed by UniProt accession numbers.
     * @param session
     * @return
     * @throws Exception
     */
    private Map<String, Protein> loadProteinsInDB(Session session) throws Exception {
        String queryString = "FROM " + Protein.class.getName();
        Query query = session.createQuery(queryString);
        List list = query.list();
        Map<String, Protein> accToProtein = new HashMap<String, Protein>();
        for (Iterator it = list.iterator(); it.hasNext();) {
            Protein protein = (Protein) it.next();
            accToProtein.put(protein.getPrimaryAccession(), 
                             protein);
        }
        return accToProtein;
    }
    
    /**
     * This method will save predicated FIs into the database based on the hibernate API. 
     * Note: 
     * 1). Only predicted FIs with a certain cutoff values have been stored in the database.
     * Not all of pairs used for predictions.
     * 2). Make sure all used protein accession numbers are from UniProt, from human, primary
     * accession numbers, have been normalized based on AA sequences. 
     * @throws Exception
     */
    @Test
    public void dumpPredictedFIs() throws Exception {
        RFPredictionResultAnalyzer rfHelper = new RFPredictionResultAnalyzer();
        // As of 2/6/2023, switch to use predicted FIs from the random forest model
        Map<String, Evidence> predictedFI2Evidence = rfHelper.loadPredictedFIsWithEvidence();
        logger.info("Total predicted FIs with evidence: " + predictedFI2Evidence.size());
        // NB by G.W (2/7/2023): One gene may be mapped to more than one UniProt accessions. However, the utilities
        // of ReactomeFIViz and the FI network are for genes. Therefore, only one UniProt accession
        // is used in this mapping.
        Map<String, String> fisInGenes2UniProt = rfHelper.mapFIsInGenesToFIsInUniProt(predictedFI2Evidence.keySet());
        // Check if any proteins cannot be mapped
        initSession();
        Session session = sessionFactory.openSession();
        // We need transaction since we need to save them
        Transaction tx = session.beginTransaction();
        // Check how many proteins should be loaded.
        // Since all UniProt accession numbers should have been sequence normalized. We can use
        // accession numbers directly.
        Map<String, Protein> accToProteinInDb = loadProteinsInDB(session);
        logger.info("Total protein from db: " + accToProteinInDb.size());
        // Help to get Proteins if not in the database
        ReactomeFuncInteractionExtractor fiHelper = new ReactomeFuncInteractionExtractor();
        List<Interaction> predictedFIs = new ArrayList<Interaction>();
        // Add these pairs to the database
        for (String pair : predictedFI2Evidence.keySet()) {
            // Switch from genes to uniprots
            String fiInUniProt = fisInGenes2UniProt.get(pair);
            int index = fiInUniProt.indexOf("\t");
            String acc1 = fiInUniProt.substring(0, index);
            String acc2 = fiInUniProt.substring(index + 1);
            Protein protein1 = accToProteinInDb.get(acc1);
            if (protein1 == null)
                protein1 = fiHelper.getProtein("UniProt:" + acc1); // All accessions are from UniProt
            if (protein1 == null)
                throw new IllegalStateException(acc1 + " cannot be found in UniProt!");
            Protein protein2 = accToProteinInDb.get(acc2);
            if (protein2 == null)
                protein2 = fiHelper.getProtein("UniProt:" + acc2);
            if (protein1 == null || protein2 == null)
                throw new IllegalStateException(acc2 + " cannot be found in UniProt!");
            Interaction interaction = new Interaction();
            interaction.setFirstProtein(protein1);
            interaction.setSecondProtein(protein2);
            interaction.setEvidence(predictedFI2Evidence.get(pair));
            predictedFIs.add(interaction);
        }
        logger.info("Total FIs will be added to the FI database: " + predictedFIs.size());
        // After all these steps, we can save them to the database.
        // Three steps are used to persist these new Interaction
        // First for Proteins
        int totalProtein = 0;
        for (Interaction interaction : predictedFIs) {
            Protein protein1 = interaction.getFirstProtein();
            if (protein1.getDbId() == 0) {// This is a new protein
                session.persist(protein1);
                totalProtein ++;
            }
            Protein protein2 = interaction.getSecondProtein();
            if (protein2.getDbId() == 0) {
                session.persist(protein2);
                totalProtein ++;
            }
        }
        logger.info("Save proteins: " + totalProtein);
        // Second for Evidence
        for (Interaction interaction : predictedFIs) {
            Evidence evidence = interaction.getEvidence();
            session.persist(evidence);
        }
        logger.info("Save evidences: " + predictedFIs.size());
        // Third for Interaction
        // I want to have JDBC control to increase the performance
        for (Interaction interaction : predictedFIs) {
            session.persist(interaction);
        }
        try  {
            session.flush();
            session.clear();
            tx.commit();
        }
        catch(Exception e) {
            tx.rollback();
        }
        // Don't forget to close the session.
        session.close();
    }
    
    /**
     * Add evidences to a list of Interactions.
     * @param interactions
     * @throws Exception
     */
    public void addEvidences(List<Interaction> interactions) throws Exception {
        // Initialize functional interactions
        FunctionalInteractionAnalyzer analyzer = new FunctionalInteractionAnalyzer();
        analyzer.setUp();
        // Map from interactions to pair
        Map<String, Interaction> pairToIntMap = new HashMap<String, Interaction>();
        Map<String, Value> valueMap = new HashMap<String, Value>();
        Protein firstProtein, secondProtein;
        String acc1, acc2;
        String pair;
        Value value;
        for (Interaction interaction: interactions) {
            firstProtein = interaction.getFirstProtein();
            acc1 = firstProtein.getPrimaryAccession();
            secondProtein = interaction.getSecondProtein();
            acc2 = secondProtein.getPrimaryAccession();
            pair = acc1 + " " + acc2;
            value = new Value();
            pairToIntMap.put(pair, interaction);
            valueMap.put(pair, value);
        }
        WEKADataAnalyzer wekaAnalyzer = new WEKADataAnalyzer();
        wekaAnalyzer.generateDataSet(valueMap);
        Instances dataset = new WEKAResultAnalyzer().createDataSet();
        double[] probs ;
        for (Iterator<String> it = valueMap.keySet().iterator(); it.hasNext();) {
            pair = it.next();
            value = valueMap.get(pair);
            probs = analyzer.calculateProbabilities(value, dataset);
            addEvidence(pairToIntMap.get(pair),
                        value, probs);
        }
    }
    
    private void addEvidence(Interaction interaction,
                             String pair,
                             Map<String, PositiveChecker> featureToChecker,
                             NaiveBayesClassifier nbc) throws Exception {
        double score = nbc.calculateScore(pair, featureToChecker);
        Evidence evidence = new Evidence();
        evidence.setScore(score);
        // For other information
        for (String feature : featureToChecker.keySet()) {
            PositiveChecker checker = featureToChecker.get(feature);
            Boolean value = checker.isPositive(pair);
            String methodName = "set" + feature.substring(0, 1).toUpperCase() + feature.substring(1);
            Method setMethod = Evidence.class.getMethod(methodName, Boolean.class);
            setMethod.invoke(evidence, value);
        }
        interaction.setEvidence(evidence);
    }
    
    private void addEvidence(Interaction interaction,
                             Value value,
                             double[] probs) {
        Evidence evidence = new Evidence();
//        // Copy values from Value to evidence
//        evidence.setProbability(probs[0]);
//        evidence.setHumanInteraction(value.humanInteraction);
//        evidence.setOrthoInteraction(value.orthoInteraction);
//        evidence.setYeastInteraction(value.yeastInteraction);
//        if (value.geneExp == null)
//            evidence.setGeneExp(GeneExpressionType.UNCORRELATED);
//        else if (value.geneExp.equals("pos"))
//            evidence.setGeneExp(GeneExpressionType.POSITIVE);
//        else if (value.geneExp.equals("neg"))
//            evidence.setGeneExp(GeneExpressionType.NEGATIVE);
//        else
//            evidence.setGeneExp(GeneExpressionType.UNCORRELATED);
//        evidence.setGoBPSemanticSimilarity(value.goBPSemSimilarity);
        interaction.setEvidence(evidence);
    }
}
