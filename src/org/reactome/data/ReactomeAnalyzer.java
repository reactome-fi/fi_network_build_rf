/*
 * Created on May 3, 2006
 *
 */
package org.reactome.data;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.log4j.Logger;
import org.gk.model.GKInstance;
import org.gk.model.InstanceUtilities;
import org.gk.model.PersistenceAdaptor;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.persistence.XMLFileAdaptor;
import org.gk.schema.InvalidAttributeException;
import org.gk.schema.SchemaAttribute;
import org.gk.schema.SchemaClass;
import org.junit.Test;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.FileUtility;
import org.reactome.fi.util.InteractionUtilities;
import org.reactome.fi.util.ReactomeUtilities;
import org.reactome.fi.util.Value;
import org.reactome.weka.FeatureHandlerForV3;

@SuppressWarnings("unchecked")
public class ReactomeAnalyzer {
	private static final Logger logger = Logger.getLogger(ReactomeAnalyzer.class);
    // Used to control if complex should be used
    protected boolean excludeComplex = false;
    protected PersistenceAdaptor dba;
    protected Long dataSourceId;
    // A helper class
    private ReactomeAnalyzerTopicHelper topicHelper;
    
    public ReactomeAnalyzer() {
        topicHelper = new ReactomeAnalyzerTopicHelper();
    }
    
    public ReactomeAnalyzerTopicHelper getTopicHelper() {
    	return this.topicHelper;
    }
    
    /**
     * Somehow there is an old gene symbol, ANP32C, for gene ANP32CP used in some FIs. Have to manually 
     * update ReferenceGeneProduct in the database (DB_ID 49306) by adding ANP32C as a gene name so that 
     * synonym mapping can work without throwing null.
     * @throws Exception
     */
    @Test
    public void fixANP32CIssueForReferenceGeneProduct() throws Exception {
        MySQLAdaptor dba = (MySQLAdaptor) getMySQLAdaptor();
        // This should be the same since this is an  old instance
        Long dbId = 49306L;
        GKInstance anp32 = dba.fetchInstance(dbId);
        if (anp32 == null) {
            logger.error("Cannot find instance for " + dbId + ".");
            return;
        }
        List<String> geneNames = anp32.getAttributeValuesList(ReactomeJavaConstants.geneName);
        String synonym = "ANP32C";
        if (geneNames.contains(synonym)) {
            logger.error(synonym + " is in the geneName slot already.");
            return;
        }
        // Make a change now
        geneNames.add(synonym);
        anp32.setAttributeValue(ReactomeJavaConstants.geneName, geneNames);
        // Make sure just update this specific attribute. Otherwise, we may have big trouble.
        // Start with transaction
        try {
            logger.info("Update the gene name of " + anp32);
            dba.startTransaction();
            dba.updateInstanceAttribute(anp32, ReactomeJavaConstants.geneName);
            dba.commit();
            // Note: We have not added IE for this update!
            logger.info("Done");
        }
        catch(Exception e) {
            logger.error("Error: " + e.getMessage(), e);
            dba.rollback();
        }
    }
    
    /**
     * Make sure dataSourceIds are correct for the version the database.
     * @return
     */
    public static List<ReactomeAnalyzer> getPathwayDBAnalyzers() {
        List<ReactomeAnalyzer> analyzers = new ArrayList<ReactomeAnalyzer>();
        // Always include ReactomeAnalyzer
        ReactomeAnalyzer reactomeAnalyzer = new ReactomeAnalyzer();
        PersistenceAdaptor dba = null;
        try {
            dba = reactomeAnalyzer.getMySQLAdaptor();
        }
        catch(Exception e) {
            System.err.println(e);
        }
        if (dba == null)
            throw new IllegalStateException("Cannot connect to a database!");
        analyzers.add(reactomeAnalyzer);
        PantherAnalyzer pantherAnalyzer = new PantherAnalyzer();
        Long pantherId = fetchDatasourceId("pantherdb", null, dba);
        pantherAnalyzer.setDataSourceId(pantherId);
        analyzers.add(pantherAnalyzer);
        // INOH is not used in version 3.
        //analyzers.add(new INOHAnalyzer());
        List<Long> dataSourceIds = new ArrayList<Long>();
        Long sourceId = fetchDatasourceId("Pathway Interaction Database", null, dba);
        ReactomeAnalyzer analyzer = new AnnotatedNCIPICAnalyzer();
        analyzer.setDataSourceId(sourceId);
        analyzers.add(analyzer);
//        dataSourceIds.add(sourceId);
        sourceId = fetchDatasourceId("BioCarta - Imported by PID", null, dba);
        dataSourceIds.add(sourceId);
        // No KEGG as of FI_2024
//        sourceId = fetchDatasourceId(KeggToReactomeConverter.KEGG_PATHWAY_DB_NAME,
//                                     KeggToReactomeConverter.KEGG_PATHWAY_URL,
//                                     dba);
//        dataSourceIds.add(sourceId);
        for (Long dataSourceId : dataSourceIds) {
            ReactomeAnalyzer tmp = new CPathAnalyzer();
            tmp.setDataSourceId(dataSourceId);
            analyzers.add(tmp);
        }
        // Add targeted interactions (TF/Target from TRED)
        TargetedInteractionAnalyzer tredAnalyzer = new TargetedInteractionAnalyzer();
        Long tredId = fetchDatasourceId("TRED", null, dba);
        tredAnalyzer.setDataSourceId(tredId);
        analyzers.add(tredAnalyzer);

        // Used a customized TargetedInteractionAnalyzer so that interactions in 
        // from ENCODE can be filtered.
        TargetedInteractionAnalyzer encodeAnalyzer = new EncodeInteractionAnalyzer();
        Long encodeId = fetchDatasourceId("ENCODE", null, dba);
        encodeAnalyzer.setDataSourceId(encodeId);
        analyzers.add(encodeAnalyzer);
        return analyzers;
    }
    
    private static Long fetchDatasourceId(String dbName,
                                          String url,
                                          PersistenceAdaptor dba) {
        Long dbId = null;
        try {
            Collection<GKInstance> instances = dba.fetchInstanceByAttribute(ReactomeJavaConstants.ReferenceDatabase,
                                                                            ReactomeJavaConstants._displayName,
                                                                            "=",
                                                                            dbName);
            if (instances != null) {
                if (url == null && instances.size() == 1) {
                    dbId = instances.iterator().next().getDBID();
                }
                else if (url != null) {
                    // In case more than one RerenceDatabases imported from different sources,
                    // since they are not normalized.
                    for (GKInstance inst : instances) {
                        String url1 = (String) inst.getAttributeValue(ReactomeJavaConstants.url);
                        if (url.equals(url1)) {
                            dbId = inst.getDBID();
                            break;
                        }
                            
                    }
                }
            }
        }
        catch(Exception e) {
            System.err.println(e);
        }
        if (dbId == null)
            throw new IllegalStateException("Cannot find a ReferenceDatabase for \"" + dbName + "\"!");
        return dbId;
    }
    
    public static String getSourceLetter(ReactomeAnalyzer analyzer) throws Exception {
        GKInstance dataSource = analyzer.getDataSource();
        return InteractionUtilities.getPathwayDBSourceLetter(dataSource);
    }
    
    public void setExcludeComplex(boolean value) {
        this.excludeComplex = value;
    }
    
    public void setMySQLAdaptor(MySQLAdaptor dba) {
        this.dba = dba;
    }
    
    public PersistenceAdaptor getMySQLAdaptor() throws Exception {
        if (dba == null) {
//            dba = new MySQLAdaptor("localhost",
//                                   "gk_central_101606",//"panther_from_david", 
//                                   //"test_reactome_plus_i_v2",
//                                   "root",
//                                   "macmysql01",
//                                   3306);
//            dba = new MySQLAdaptor("localhost",
//                                   "reactome_28_plus_i",//"panther_from_david", 
//                                   "root",
//                                   "macmysql01",
//                                   3306);
            dba = new MySQLAdaptor("localhost",
                                   FIConfiguration.getConfiguration().get("REACTOME_SOURCE_DB_NAME"), 
                                   FIConfiguration.getConfiguration().get("DB_USER"),
                                   FIConfiguration.getConfiguration().get("DB_PWD"),
                                   3306);
//            dba = new MySQLAdaptor("localhost",
//                                   "gk_central_031309",//"panther_from_david", 
//                                   "root",
//                                   "macmysql01",
//                                   3306);
        }
        return dba;
    }
    
    public void setDataSourceId(Long id) {
        this.dataSourceId = id;
    }
    
    public Set<String> generateUniProtPairsFromTopics() throws Exception {
        Set<String> pairs = new HashSet<String>();
        Map<GKInstance, Set<String>> topics2Ids = grepIDsFromTopics();
        for (Iterator<GKInstance> it = topics2Ids.keySet().iterator(); it.hasNext();) {
            GKInstance topic = it.next();
            Set<String> ids = topics2Ids.get(topic);
            List<String> idList = new ArrayList<String>(ids);
            Collections.sort(idList);
            for (int i = 0; i < idList.size() - 1; i++) {
                String id1 = idList.get(i);
                for (int j = i + 1; j < idList.size(); j++) {
                    String id2 = idList.get(j);
                    pairs.add(id1 + " " + id2);
                }
            }
        }
        System.out.println("Total Pairs: " + pairs.size());
        return pairs;
    }
    
    /**
     * Using this method to generate interaction list for each topic.
     * @return
     * @throws Exception
     */
    public Map<GKInstance, Set<String>> grepInteractionsForTopics() throws Exception {
        Map<GKInstance, Set<String>> topicToInteraction = new HashMap<GKInstance, Set<String>>();
        List<GKInstance> topics = getTopics();
        prepareReactions();
        prepareComplexes();
        Set<GKInstance> interactors = new HashSet<GKInstance>();
        for (GKInstance topic : topics) {
            Set<String> interactions = grepInteractionsForTopic(topic);
            topicToInteraction.put(topic, interactions);
        }
        return topicToInteraction;
    }
    
    /**
     * Use this method to load a list of FIs for a specified pathway (here is used as topic).
     * @param topic
     * @return
     * @throws Exception
     */
    public Set<String> grepInteractionsForTopic(GKInstance topic) throws Exception {
        Set<GKInstance> components = topicHelper.grepPathwayComponents(topic);
        Set<String> interactions = new HashSet<String>();
        Set<GKInstance> interactors = new HashSet<GKInstance>();
        for (GKInstance comp : components) {
            if (comp.getSchemClass().isa(ReactomeJavaConstants.Reaction)) {
                interactors.clear();
                extractInteractorsFromReaction(comp, interactors);
                generateInteractions(interactors, interactions, comp);
            }
            else if (comp.getSchemClass().isa(ReactomeJavaConstants.Interaction)) {
                List list = comp.getAttributeValuesList(ReactomeJavaConstants.interactor);
                if (list != null) {
                    interactors.clear();
                    interactors.addAll(list);
                    generateInteractions(interactors, interactions, comp);
                }
            }
            else if (comp.getSchemClass().isa(ReactomeJavaConstants.Complex)) {
                interactors.clear();
                grepComplexComponents(comp, interactors);
                generateInteractions(interactors, interactions, comp);
            }
        }
        return interactions;
    }
    
    public Set<String> generateUniProtIdsFromTopics() throws Exception {
        Map<GKInstance, Set<String>> topics2Ids = grepIDsFromTopics();
        Set<String> uniProtIds = new HashSet<String>();
        for (Iterator<GKInstance> it = topics2Ids.keySet().iterator(); it.hasNext();) {
            GKInstance topic = it.next();
            Set<String> ids = topics2Ids.get(topic);
            uniProtIds.addAll(ids);
        }
        return uniProtIds;
    }
    
    protected List<GKInstance> getTopics() throws Exception {
        MySQLAdaptor releasedDBA = (MySQLAdaptor) getMySQLAdaptor();
        // The following code should be used for generating NBC training data set (May 11, 2009)
//        Collection frontPages = releasedDBA.fetchInstancesByClass(ReactomeJavaConstants.FrontPage);
//        List<GKInstance> topics = new ArrayList<GKInstance>();
//        for (Iterator it = frontPages.iterator(); it.hasNext();) {
//            GKInstance frontPage = (GKInstance) it.next();
//            List items = frontPage.getAttributeValuesList(ReactomeJavaConstants.frontPageItem);
//            for (Iterator it1 = items.iterator(); it1.hasNext();)
//                topics.add((GKInstance)it1.next());
//        }
        // As of April 19, 2007, a list of Reactome pathways is fetched from a semi-manually
        // create file, ReactomePathways.txt.
        List<GKInstance> topics = new ArrayList<GKInstance>();
        String fileName = FIConfiguration.getConfiguration().get("REACTOME_PATHWAYS");
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = null;
        int index = 0;
        Long id = null;
        while ((line = fu.readLine()) != null) {
            index = line.indexOf("\t");
            id = Long.parseLong(line.substring(0, index));
            GKInstance pathway = releasedDBA.fetchInstance(id);
            topics.add(pathway);
        }
        return topics;
    }
    
    public Map<GKInstance, Map<String, Integer>> grepIDNumberFromTopics() throws Exception {
        List<GKInstance> topics = getTopics();
        Map<GKInstance, Map<String, Integer>> topicToIDNumber = new HashMap<GKInstance, Map<String, Integer>>();
        for (GKInstance topic : topics) {
            Map<String, Integer> id2Number = topicHelper.grepIDToNumberFromTopic(topic);
            topicToIDNumber.put(topic, id2Number);
        }
        return topicToIDNumber;
    }
    
    public Map<GKInstance, Set<String>> grepIDsFromTopics() throws Exception {
        List<GKInstance> topics = getTopics();
        return grepIDsFromTopics(topics);
    }

    public Map<GKInstance, Set<String>> grepIDsFromTopics(Collection<GKInstance> topics) throws Exception {
        // Try to get ids for each topics
        long time1 = System.currentTimeMillis();
        Map<GKInstance, Set<String>> topics2Ids = new HashMap<GKInstance, Set<String>>();
        for (GKInstance topic : topics) {
//            System.out.println("Topic: " + topic);
            Set<String> ids = grepIDsFromTopic(topic);
            topics2Ids.put(topic, ids);
        }
        long time2 = System.currentTimeMillis();
        logger.info("Time for grepping IDs for topics: " + (time2 - time1));
        // Print out the id numbers in each topic.
        for (Iterator<GKInstance> it = topics2Ids.keySet().iterator(); it.hasNext();) {
            GKInstance topic = it.next();
            Set<String> ids = topics2Ids.get(topic);
            logger.info(topic.getDisplayName() + ": " + ids.size());
        }
        return topics2Ids;
    }
    
    protected Set<GKInstance> grepPathwayParticipants(GKInstance pathway) throws Exception {
        return topicHelper.grepPathwayParticipants(pathway);
    }
    
    public Set<String> grepIDsFromTopic(GKInstance topic) throws Exception {
        Set<String> ids = new HashSet<String>();
        // First load all PhysicalEntities involved in Reactions
        Set<GKInstance> participants = grepPathwayParticipants(topic);
        // Grep ReferencePeptideSequence
        for (GKInstance participant : participants) {
            if (excludeComplex && participant.getSchemClass().isa(ReactomeJavaConstants.Complex))
                continue;
            Set<GKInstance> refPeptides = grepRefPepSeqs(participant);
            for (GKInstance ref : refPeptides) {
                GKInstance refDb = (GKInstance) ref.getAttributeValue(ReactomeJavaConstants.referenceDatabase);
                // It should be faster by comparing DB_ID
                if (refDb == null || !refDb.getDBID().equals(2L))
                    continue;
                String identifier = (String) ref.getAttributeValue(ReactomeJavaConstants.identifier);
                if (identifier != null)
                    ids.add(identifier);
            }
        }
        return ids;
    }
    
    private Collection loadReactions() throws Exception {
        MySQLAdaptor dba = (MySQLAdaptor) getMySQLAdaptor();
        // Get homo sapiens
        GKInstance homosapiens = dba.fetchInstance(48887L);
        // Load all reactions for analyzed: human reactions only
        Collection reactions = dba.fetchInstanceByAttribute(ReactomeJavaConstants.Reaction,
                                                            ReactomeJavaConstants.species,
                                                            "=",
                                                            homosapiens);
        // Load precedingEvent values
        SchemaClass schema = dba.getSchema().getClassByName(ReactomeJavaConstants.Reaction);
        SchemaAttribute att = schema.getAttribute(ReactomeJavaConstants.precedingEvent);
        dba.loadInstanceAttributeValues(reactions, att);
        return reactions;
    }
    
    @SuppressWarnings("rawtypes")
    public void extractInteractorsFromReaction(GKInstance rxn, 
                                               Set<GKInstance> interactors) throws Exception {
        List input = rxn.getAttributeValuesList(ReactomeJavaConstants.input);
//        // Something special for gene regulatory reaction annotated in BlackBoxEvent
//        boolean isGeneRegulatory = false;
//        if (rxn.getSchemClass().isa(ReactomeJavaConstants.BlackBoxEvent) && input != null && input.size() == 1) {
//            GKInstance input1 = (GKInstance) input.get(0);
//            if (input1.getSchemClass().isValidAttribute(ReactomeJavaConstants.referenceEntity)) {
//                GKInstance refEntity = (GKInstance) input1.getAttributeValue(ReactomeJavaConstants.referenceEntity);
//                if (refEntity != null && refEntity.getSchemClass().isa(ReactomeJavaConstants.ReferenceDNASequence)) {
//                    isGeneRegulatory = true;
//                    GKInstance output = (GKInstance) rxn.getAttributeValue(ReactomeJavaConstants.output);
//                    if (output != null)
//                        interactors.add(output);
//                }
//            }
//        }
//        if (!isGeneRegulatory && input != null) {
//            interactors.addAll(input);
//        }
        if (input != null)
            interactors.addAll(input);
        // Get catalyst
        List cas = rxn.getAttributeValuesList(ReactomeJavaConstants.catalystActivity);
        if (cas != null) {
            for (Iterator it1 = cas.iterator(); it1.hasNext();) {
                GKInstance ca = (GKInstance) it1.next();
                List catalysts = ca.getAttributeValuesList(ReactomeJavaConstants.physicalEntity);
                if (catalysts != null)
                    interactors.addAll(catalysts);
            }
        }
        // Check regulators
        Collection regulations = InstanceUtilities.getRegulations(rxn);
        if (regulations != null) {
            for (Iterator it1 = regulations.iterator(); it1.hasNext();) {
                GKInstance regulation = (GKInstance) it1.next();
                List regulators = regulation.getAttributeValuesList(ReactomeJavaConstants.regulator);
                for (Iterator it2 = regulators.iterator(); it2.hasNext();) {
                    GKInstance regulator = (GKInstance) it2.next();
                    if (regulator.getSchemClass().isa(ReactomeJavaConstants.PhysicalEntity))
                        interactors.add(regulator);
                }
            }
        }
    }
    
    @Test
    public void extractInteractionsInGenes() throws Exception {
        Collection reactions = prepareReactions();
        Collection complexes = prepareComplexes();
        Set<String> interactions = new HashSet<String>();
        GKInstance rxn = null;
        Set<GKInstance> interactors = new HashSet<GKInstance>();
        long time1 = System.currentTimeMillis();
        for (Iterator it = reactions.iterator(); it.hasNext();) {
            rxn = (GKInstance) it.next();
            //System.out.println("Reaction: " + c++);
            extractInteractorsFromReaction(rxn, interactors);
            generateInteractions(interactors, interactions, rxn, true);
            interactors.clear();
        }
        System.out.println("Total interactions from reactions: " + interactions.size());
        GKInstance complex = null;
        for (Iterator it = complexes.iterator(); it.hasNext();) {
            complex = (GKInstance) it.next();
            //System.out.println("Complex: " + c++ + " " + complex.getDBID());
            interactors.clear();
            grepComplexComponents(complex, interactors);
            // No need
            //if (interactors.size() > 10)
            //    continue; // cutoff set manually
            generateInteractions(interactors, interactions, complex, true);
        }
        long time2 = System.currentTimeMillis();
        System.out.println("Time for looping: " + (time2 - time1));
        System.out.println("Total interactions from Reactome: " + interactions.size());
        // Perform some cleaning up
        Set<String> cleaned = new HashSet<>();
        int index = 0;
        for (String fi : interactions) {
            index = fi.indexOf("\t");
            fi = InteractionUtilities.generateFIFromGene(fi.substring(0, index),
                                                         fi.substring(index + 1));
            cleaned.add(fi);
        }
        System.out.println("Total FIs after reordering: " + cleaned.size());
        FileUtility fu = new FileUtility();
        String fileName = FIConfiguration.getConfiguration().get("REACTOME_FI_FILE");
        index = fileName.lastIndexOf(".");
        fileName = fileName.substring(0, index) + "_Genes.txt";
        fu.saveInteractions(cleaned, fileName);
    }
    
    /**
     * This method is used to dump FIs extracted from the Reactome database to 
     * an external file.
     * @throws Exception
     */
    @Test
    public void extractInteractions() throws Exception {
        Set<String> interactionSet = extractInteractionSet();
        FileUtility fu = new FileUtility();
        //String fileName = "results/v2/ReactomeInteractions020507.txt";
        // 28 for the release number
        // Do a filters
        ProteinIdFilters filters = new ProteinIdFilters();
        System.out.println("Before filtering: " + interactionSet.size());
        Set<String> filtered = filters.cleanUpVsUniProt(interactionSet);
        System.out.println("After filtering: " + filtered.size());
//        String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR + "ReactomeInteractions28_051711.txt";
//        String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR + "ReactomeInteractions36.txt";
        String fileName = FIConfiguration.getConfiguration().get("REACTOME_FI_FILE");
        fu.saveInteractions(filtered, fileName);
        // Check what have been removed
        // Note: interactions filtered out are extracted from reactions involed other non-human
        // species (e.g. virus, bacteria, etc).
//        interactionSet.removeAll(filtered);
//        int count = 0;
//        System.out.println("Total removed: " + interactionSet.size());
//        for (String interaction : interactionSet) {
//            System.out.println(interaction);
//            count ++;
//            if (count == 10)
//                break;
//        }
    }
    
    @Test
    public void checkNumberChangesBetweenTwoReleases() throws Exception {
        String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR") + "ReactomeOnlyInteractions28_051711.txt";
        FileUtility fu = new FileUtility();
        Set<String> fis1 = fu.loadInteractions(fileName);
        Set<String> fiProteins1 = InteractionUtilities.grepIDsFromInteractions(fis1);
        fileName = FIConfiguration.getConfiguration().get("RESULT_DIR") + "ReactomeOnlyInteractions36.txt";
        Set<String> fis2 = fu.loadInteractions(fileName);
        Set<String> fiProteins2 = InteractionUtilities.grepIDsFromInteractions(fis2);
        System.out.println("Release 28 FIs: " + fis1.size());
        System.out.println("           Proteins: " + fiProteins1.size());
        System.out.println("Release 36 FIs: " + fis2.size());
        System.out.println("           Proteins: " + fiProteins2.size());
        // Check how many new reactions have been predicted before
        Set<String> totalFIs = fu.loadInteractions(FIConfiguration.getConfiguration().get("INTERACTION_FILE_NAME"));
        Set<String> totalFIProteins = InteractionUtilities.grepIDsFromInteractions(totalFIs);
        System.out.println("Total FIs: " + totalFIs.size());
        Set<String> newFIs = new HashSet<String>(fis2);
        newFIs.removeAll(fis1);
        Set<String> shared = InteractionUtilities.getShared(totalFIs, newFIs);
        Set<String> sharedProteins = InteractionUtilities.grepIDsFromInteractions(shared);
        System.out.println("Total new FIs: " + newFIs.size());
        System.out.println("    Shared FIs: " + shared.size());
        System.out.println("    Proteins: " + sharedProteins.size());
        Set<String> newProteins = new HashSet<String>(fiProteins2);
        newProteins.removeAll(fiProteins1);
        shared = InteractionUtilities.getShared(totalFIProteins, newProteins);
        System.out.println("Total new proteins: " + newProteins.size());
        System.out.println("Shared new Proteins: " + shared.size());
//        MySQLAdaptor dba = (MySQLAdaptor) getMySQLAdaptor();
//        Collection<?> c = dba.fetchInstanceByAttribute(ReactomeJavaConstants.ReferenceGeneProduct, 
//                                                                     ReactomeJavaConstants.species, 
//                                                                     "=",
//                                                                     48887L);
//        dba.loadInstanceAttributeValues(c, new String[]{ReactomeJavaConstants.identifier});
//        Set<String> idsInDb = new HashSet<String>();
//        for (Object obj : c) {
//            GKInstance inst = (GKInstance) obj;
//            String id = (String) inst.getAttributeValue(ReactomeJavaConstants.identifier);
//            idsInDb.add(id);
//        }
//        System.out.println("Total ids in db: " + idsInDb.size());
//        Set<String> notInFIs = new HashSet<String>(idsInDb);
//        notInFIs.removeAll(fiProteins2);
//        System.out.println("Not in FIs: " + notInFIs.size());
//        int count = 0;
//        for (String id : notInFIs) {
//            System.out.println(id);
//            count ++;
//            if (count > 100)
//                break;
//        }
    }
    
    /**
     * Use this method to load a list of pre-generated FIs from the Reactome database.
     * @return
     * @throws IOException
     */
    public Set<String> loadFIsFromFile() throws IOException {
        FileUtility fu = new FileUtility();
        //String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR + "ReactomeInteractions28.txt";
//        String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR + "FIs_Reactome.txt";
        String fileName = FIConfiguration.getConfiguration().get("REACTOME_FI_FILE");
        return fu.loadInteractions(fileName);
    }
    
    public Set<String> extractInteractionSetWithComplexAndSet() throws Exception {
        // Extract from gk_central reactions
        Collection reactions = prepareReactions();
        Collection complexes = prepareComplexes();
        Set<String> interactions = new HashSet<String>();
        GKInstance rxn = null;
        Set<GKInstance> interactors = new HashSet<GKInstance>();
        long time1 = System.currentTimeMillis();
        int c = 0;
        for (Iterator it = reactions.iterator(); it.hasNext();) {
            rxn = (GKInstance) it.next();
            //System.out.println("Reaction: " + c++);
            extractInteractorsFromReaction(rxn, interactors);
            generateInteractions(interactors, interactions, rxn);
            interactors.clear();
        }
        System.out.println("Total interactions from " + getDataSource() + ": "
                           + interactions.size());
        return interactions;
    }
    
    protected void generateInteractionsWithComplexAndSet(Set<GKInstance> interactors,
                                                         Set<String> interactions,
                                                         GKInstance reaction) throws Exception {
        List<GKInstance> list = new ArrayList<GKInstance>(interactors);
        int size = list.size();
        for (int i = 0; i < size - 1; i++) {
            GKInstance interactor1 = list.get(i);
            for (int j = i + 1; j < size; j++) {
                GKInstance interactor2 = list.get(j);
                generateInteractionsWithComplexAndSet(interactor1, 
                                                      interactor2, 
                                                      interactions);
            }
        }
    }
    
    private void generateInteractionsWithComplexAndSet(GKInstance interactor1,
                                                       GKInstance interactor2,
                                                       Set<String> interactions) throws Exception {
        Set<GKInstance> refPepSeqs1 = grepRefPepSeqs(interactor1);
        if (refPepSeqs1.size() == 0)
            return;
        Set<GKInstance> refPepSeqs2 = grepRefPepSeqs(interactor2);
        if (refPepSeqs2.size() == 0)
            return;
        String int1 = convertRefPepSeqToString(refPepSeqs1, interactor1);
        if (int1.length() == 0)
            return;
        String int2 = convertRefPepSeqToString(refPepSeqs2, interactor2);
        if (int2.length() == 0)
            return;
        int compare = int1.compareTo(int2);
        if (compare < 0)
            interactions.add(int1 + " " + int2);
        else
            interactions.add(int2 + " " + int1);
    }
              
    private String convertRefPepSeqToString(Set<GKInstance> refPepSeqs,
                                            GKInstance interactor) throws Exception {
        List<String> ids = new ArrayList<String>();
        for (Iterator<GKInstance> it = refPepSeqs.iterator(); it.hasNext();) {
            GKInstance refPepSeq = it.next();
            String identifier = (String) refPepSeq.getAttributeValue(ReactomeJavaConstants.identifier);
            if (identifier == null)
                continue; // maybe
            ids.add(identifier);
        }
        Collections.sort(ids);
        String delimit = "?"; // default: should not be used
        if (interactor.getSchemClass().isa(ReactomeJavaConstants.Complex))
            delimit = ",";
        else if (interactor.getSchemClass().isa(ReactomeJavaConstants.EntitySet))
            delimit = "|";
        StringBuilder builder = new StringBuilder();
        for (Iterator<String> it = ids.iterator(); it.hasNext();) {
            builder.append(it.next());
            if (it.hasNext())
                builder.append(delimit);
        }
        return builder.toString();
    }
    
    private Set<GKInstance> grepComplexes(Set<GKInstance> interactors) throws Exception {
        Set<GKInstance> complexes = new HashSet<GKInstance>();
        for (GKInstance interactor : interactors) {
            if (interactor.getSchemClass().isa(ReactomeJavaConstants.Complex))
                complexes.add(interactor);
            Set<GKInstance> temp = InstanceUtilities.getContainedInstances(interactor,
                                                                           ReactomeJavaConstants.hasMember,
                                                                           ReactomeJavaConstants.hasComponent);
            for (GKInstance inst : temp) {
                if (inst.getSchemClass().isa(ReactomeJavaConstants.Complex))
                    complexes.add(inst);
            }
        }
        return complexes;
    }
    
    /**
     * Generate a mapping file from genes to reactions.
     * @throws Exception
     */
    @Test
    public void generateGenesToReactionsMap() throws Exception {
        Collection reactions = prepareReactions();
        Collection complexes = prepareComplexes();
        // Get the directory
        String fileName = FIConfiguration.getConfiguration().get("GENE_TO_REACTOME_REACTIONS");
//        String fileName = resultDir + File.separator + "ReactomeGenesToReactions_082316.txt";
        FileUtility fu = new FileUtility();
        fu.setOutput(fileName);
        // No header is required
        for (Iterator it = reactions.iterator(); it.hasNext();) {
            GKInstance rxn = (GKInstance) it.next();
            Set<String> genes = grepGenesFromReaction(rxn);
            if (genes.size() == 0)
                continue;
            for (String gene : genes)
                fu.printLine(gene + "\t" + rxn.getDisplayName());
        }
        fu.close();
    }
    
    public Set<String> grepGenesFromComplex(GKInstance complex) throws Exception {
        Set<String> genes = new HashSet<String>();
        ReactomeUtilities.grepGenesFromEntity(complex, genes);
        return genes;
    }
    
    public Set<String> grepGenesFromReaction(GKInstance rxn) throws Exception {
        Set<GKInstance> participants = InstanceUtilities.getReactionParticipants(rxn);
        Set<String> genes = new HashSet<String>();
        for (GKInstance participant : participants) {
            ReactomeUtilities.grepGenesFromEntity(participant, genes);
        }
        return genes;
    }
    
    @Test
    public void generateFIsForOneReaction() throws Exception {
        MySQLAdaptor dba = (MySQLAdaptor) getMySQLAdaptor();
        Long dbId = 5672965L; // RAS GEFs promote RAS nucleotide exchange
//        dbId = 5672972L; // MAP2Ks and MAPKs bind to the activated RAF complex
//        dbId = 69213L; // Formation of Cyclin D:Cdk4/6 complexes 
//        dbId = 5617896L; // Retinoic acid activates HOXD4 chromatin
        String dirName = "/Users/gwu/Documents/EclipseWorkspace/caBigR3/results/DriverGenes/Drivers_0816/";
        String fileName = dirName + "FIsInReaction" + dbId + ".txt";
        generateFIsForOneReaction(dba,
                                  dbId,
                                  fileName);
    }
    
    /**
     * Extract Reactome FIs for complexes and reactions related to cancer genes. 
     * @throws Exception
     */
    @Test
    public void generateFIsForCancerReactionsAndComplexes() throws Exception {
        String dirName = "/Users/gwu/Documents/EclipseWorkspace/caBigR3/results/DriverGenes/Drivers_0816/";
        // Load known missense cancer driver genes
        Set<String> cancerDrivers = new HashSet<String>();
        String cancerGeneFileName = dirName + "AllKnownMissenseDriver.txt";
        FileUtility fu = new FileUtility();
        fu.setInput(cancerGeneFileName);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            cancerDrivers.add(tokens[0]);
        }
        fu.close();
        System.out.println("Total cancer drivers: " + cancerDrivers.size());
        
        // Grep FIs from complexes and reactions
        Set<String> fis = new HashSet<String>();
        Collection<GKInstance> complexes = prepareComplexes();
        System.out.println("Total complexes: " + complexes.size());
        Collection<GKInstance> reactions = prepareReactions();
        System.out.println("Total reactions: " + reactions.size());
        Set<GKInstance> cancerReactions = new HashSet<GKInstance>();
        Set<GKInstance> cancerComplexes = new HashSet<GKInstance>();
        for (GKInstance reaction : reactions) {
            Set<String> genes = grepGenesFromReaction(reaction);
            if (genes == null || genes.size() == 0)
                continue;
            if (InteractionUtilities.isShared(genes, cancerDrivers))
                cancerReactions.add(reaction);
        }
        for (GKInstance complex : complexes) {
            Set<String> genes = grepGenesFromComplex(complex);
            if (genes == null || genes.size() == 0)
                continue;
            if (InteractionUtilities.isShared(genes, cancerDrivers))
                cancerComplexes.add(complex);
        }
        System.out.println("Total cancer reactions: " + cancerReactions.size());
        System.out.println("Total cancer complexes: " + cancerComplexes.size());
        Set<String> cancerFIs = extractInteractionSet(cancerReactions, cancerComplexes);

        // Export these FIs
        String outFileName = dirName + "FIsInCancerReactionsAndComplexes_122116.txt";
        fu.saveCollection(cancerFIs, outFileName);
    }
    
    @Test
    public void generateFIsForReactionsWithFeatures() throws Exception {
        MySQLAdaptor dba = (MySQLAdaptor) getMySQLAdaptor();
//        String dirName = "/Users/gwu/Documents/EclipseWorkspace/caBigR3/results/DriverGenes/Drivers_0816/";
//        
//        String rxtFile = dirName + "SelectedForInteractome3d_091416.txt";
//        List<Long> dbIds = loadReactionIds(rxtFile);
//        System.out.println("Total DB_IDs: " + dbIds.size());
//        List<Long> dbIdsFDR05 = loadReactionIds(dirName + "SelectedForInteractome3d_fdr_05_091416.txt");
//        System.out.println("Total DB_IDs for FDR <= 0.05: " + dbIdsFDR05.size());
//        
//        String outFileName = dirName + "FIsInSelectedReactions_FDR_05_092516.txt";
        
        // Get a set of DB_IDs 
        String fileName = FIConfiguration.getConfiguration().get("FI_TO_REACTOME_REACTIONS");
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        Set<Long> dbIds = new HashSet<Long>();
        String line = null;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            dbIds.add(new Long(tokens[2]));
        }
        fu.close();
        
        String outFileName = FIConfiguration.getConfiguration().get("RESULT_DIR") + File.separator + "ReactomeFIsWithPPIFeatures_041117.txt";
        generateFIsForReactionsWithFeatures(dba, 
                                            new ArrayList<Long>(dbIds),
                                            outFileName);
        
        
//        dbIdsFDR05.removeAll(dbIds); // Reactions between [0.05, 0.01)
//        System.out.println("Total DB_IDs for FDRs [0.05, 0.01): " + dbIdsFDR05.size());
//        generateFIsForReactionsWithFeatures(dba,
//                                            dbIdsFDR05, 
//                                            dirName + "FIsInSelectedForInteractome3d_FDR_05_01_091416.txt");
        
//        String outFileName = dirName + "FIsInSelectedForInteractome3d_091416.txt";
//        generateFIsForReactionsWithFeatures(dba, dbIds, outFileName);
    }

    private List<Long> loadReactionIds(String rxtFile) throws IOException {
        FileUtility fu = new FileUtility();
        fu.setInput(rxtFile);
        String line = fu.readLine();
        List<Long> dbIds = new ArrayList<Long>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            dbIds.add(new Long(tokens[0]));
        }
        fu.close();
        return dbIds;
    }
    
    private void generateFIsForOneReaction(MySQLAdaptor dba,
                                           Long dbId,
                                           String fileName) throws Exception, IOException {
        List<Long> dbIds = new ArrayList<Long>();
        dbIds.add(dbId);
        generateFIsForReactionsWithFeatures(dba, dbIds, fileName);
    }

    public void generateFIsForReactionsWithFeatures(MySQLAdaptor dba,
                                                    Collection<Long> dbIds,
                                                    String fileName) throws Exception, IOException {
        Set<String> interactions = new HashSet<String>();
        Set<String> rxtFIs = new HashSet<String>();
        Set<GKInstance> interactors = new HashSet<GKInstance>();
        Map<String, Set<Long>> fiToRxtIDs = new HashMap<String, Set<Long>>();
        for (Long dbId : dbIds) {
            GKInstance reaction = dba.fetchInstance(dbId);
            interactors.clear();
            rxtFIs.clear();
            generateFIsForSingleReaction(interactors, rxtFIs, reaction, false, false);
            if (rxtFIs.size() == 0)
                continue;
            interactions.addAll(rxtFIs);
            for (String fi : rxtFIs) {
                InteractionUtilities.addElementToSet(fiToRxtIDs,
                                                     fi, 
                                                     dbId);
            }
        }
        
        // Feature handler
        FeatureHandlerForV3 featureHandler = new FeatureHandlerForV3();
        List<String> featureList = featureHandler.getFeatureList();
        Map<String, Value> fiToValue = featureHandler.convertPairsToValues(interactions, true);
        
        Map<String, String> idToGene = getUniProtToGeneMap(dba);
        
//        for (String fi : interactions)
//            System.out.println(fi.replace('\t', ' '));
        FileUtility fu = new FileUtility();
        fu.setOutput(fileName);
        fu.printLine("UniProt1\tUniProt2\tGene1\tGene2\tPositiveFeature\tHumanPPI\tMousePPI\tFlyPPI\tWormPPI\tYeastPPI\tDomainInt\tReactions");
        for (String fi : interactions) {
            Value value = fiToValue.get(fi);
            Boolean posFeature = value.humanInteraction |
                                 value.mousePPI |
                                 value.dmePPI |
                                 value.celPPI |
                                 value.scePPI |
                                 value.pfamDomainInt;
            int index = fi.indexOf("\t");
            String gene1 = idToGene.get(fi.substring(0, index));
            String gene2 = idToGene.get(fi.substring(index + 1));
            String rxtIds = fiToRxtIDs.get(fi).toString();
            rxtIds = rxtIds.substring(1, rxtIds.length() - 1);
            fu.printLine(fi + "\t" + 
                               gene1 + "\t" + 
                               gene2 + "\t" +
                               posFeature + "\t" + 
                               value.humanInteraction + "\t" + 
                               value.mousePPI + "\t" + 
                               value.dmePPI + "\t" + 
                               value.celPPI + "\t" + 
                               value.scePPI + "\t" +
                               value.pfamDomainInt + "\t" + 
                               rxtIds);
        }
        fu.close();
    }
    
    @Test
    public void testGeneratePathwayIds() throws Exception {
    	PersistenceAdaptor dba = getMySQLAdaptor();
        // For a single test
        GKInstance pathway = dba.fetchInstance(9609736L);
        getTopicHelper().setNeedCandidateRepeatedUnit(true);
        Set<String> ids = grepIDsFromTopic(pathway);
        System.out.println("Total ids: " + ids.size());
        ids.stream().sorted().forEach(System.out::println);
    }
    
    @Test
    public void dumpReactomeHumanPathways() throws Exception {
        PersistenceAdaptor dba = getMySQLAdaptor();
        GKInstance homosapiens = dba.fetchInstance(48887L);
        Collection<GKInstance> pathways = dba.fetchInstanceByAttribute(ReactomeJavaConstants.Pathway,
                ReactomeJavaConstants.species,
                "=", 
                homosapiens);
       System.out.println("Total human pathways: " + pathways.size());
       
       String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR") + File.separator +  "ReactomeHumanPathways_04292024.txt";
       FileUtility fu = new FileUtility();
       fu.setOutput(fileName);
       fu.printLine("dbId\tDisplayName");
       for (GKInstance pathway : pathways) {
           fu.printLine(pathway.getDBID() + "\t" + pathway.getDisplayName());
       }
       fu.close();
    }
    
    @Test
    public void testGenerateFIsForSingleReaction() throws Exception {
        PersistenceAdaptor dba = getMySQLAdaptor();
        // For a single test
        GKInstance reaction = dba.fetchInstance(6802941L);
        Set<GKInstance> interactors = new HashSet<GKInstance>();
        Set<String> interactions = new HashSet<String>();
        generateFIsForSingleReaction(interactors, interactions, reaction, false, true);
        System.out.println("Total interactions: " + interactions.size());
    }
    
    @Test
    public void generateFIsToReactionsMap() throws Exception {
        Collection reactions = prepareReactions();
        Collection complexes = prepareComplexes();
        // Get the directory
        String fileName = FIConfiguration.getConfiguration().get("FI_TO_REACTOME_REACTIONS");
//        String fileName = resultDir + File.separator + "ReactomeFIsToReactions_082216.txt";
//        String fileName = resultDir + File.separator + "ReactomeFIsToReactionsWithComplexes_082516.txt";
        FileUtility fu = new FileUtility();
        fu.setOutput(fileName);
        // No header is required
        Set<GKInstance> interactors = new HashSet<GKInstance>();
        Set<String> interactions = new HashSet<String>();
        for (Iterator it = reactions.iterator(); it.hasNext();) {
            GKInstance rxn = (GKInstance) it.next();
            //System.out.println("Reaction: " + c++);
            generateFIsForSingleReaction(interactors, interactions, rxn, false, true);
            interactors.clear();
            if (interactions.size() == 0)
                continue;
            for (String interaction : interactions) {
                fu.printLine(interaction + "\t" +  rxn.getDBID() + "\t" + rxn.getDisplayName());
            }
            interactions.clear();
            interactors.clear();
        }
        fu.close();
    }

    private void generateFIsForSingleReaction(Set<GKInstance> interactors,
                                              Set<String> interactions,
                                              GKInstance rxn,
                                              boolean includeFIsInComplex,
                                              boolean useGeneName) throws Exception {
        extractInteractorsFromReaction(rxn, interactors);
        generateInteractions(interactors, interactions, rxn, useGeneName);
        if (!includeFIsInComplex)
            return;
        // Collect FIs from complexes involved in reactions
        Set<GKInstance> rxnComplexes = grepComplexes(interactors);
        for (GKInstance complex : rxnComplexes) {
            Set<GKInstance> complexInteractors = new HashSet<GKInstance>();
            grepComplexComponents(complex, complexInteractors);
            generateInteractions(complexInteractors, 
                                 interactions,
                                 complex,
                                 useGeneName);
        }
    }
    
    /**
     * Generate a map from reactions to pathways.
     * @throws Exception
     */
    @Test
    public void generateReactionsToPathwaysMap() throws Exception {
        Collection reactions = prepareReactions();
        MySQLAdaptor dba = (MySQLAdaptor) getMySQLAdaptor();
        Collection<GKInstance> pathways = dba.fetchInstanceByAttribute(ReactomeJavaConstants.Pathway,
                                                              ReactomeJavaConstants.species,
                                                              "=",
                                                              48887);
        dba.loadInstanceAttributeValues(pathways, new String[]{
                ReactomeJavaConstants.hasEvent,
                ReactomeJavaConstants.dataSource,
                ReactomeJavaConstants.disease
        });
        System.out.println("Total human pathways: " + pathways.size());
        // Remove disease pathway
        Long diseaseId = 1643685L;
        GKInstance disease = dba.fetchInstance(diseaseId);
        Set<GKInstance> diseasePathways = InstanceUtilities.getContainedEvents(disease);
        for (GKInstance diseasePathway : diseasePathways) {
            // Normal pathways may be listed under disease in the old approach for drawing pathway diagrams
            GKInstance diseaseAttribute = (GKInstance) diseasePathway.getAttributeValue(ReactomeJavaConstants.disease);
            if (diseaseAttribute != null)
                pathways.remove(diseasePathway);
        }
        pathways.remove(disease);
        System.out.println("Remove disease pathways: " + pathways.size());
        // Get the directory
        String resultDir = FIConfiguration.getConfiguration().get("RESULT_DIR");
//        String fileName = resultDir + File.separator + "ReactomeReactionsToPathways_090116.txt";
        String fileName = resultDir + File.separator + "ReactomeReactionsToPathways_051017.txt";
        FileUtility fu = new FileUtility();
        fu.setOutput(fileName);
        // No header is required
        for (GKInstance pathway : pathways) {
            if(pathway.getAttributeValue(ReactomeJavaConstants.dataSource) != null)
                continue; // Use Reactome only
            System.out.println("Handling " + pathway);
            Set<GKInstance> pathwayEvents = InstanceUtilities.grepPathwayEventComponents(pathway);
            for (GKInstance event : pathwayEvents) {
                if (!(event.getSchemClass().isa(ReactomeJavaConstants.ReactionlikeEvent)))
                    continue;
                fu.printLine(event.getDBID() + "\t" +  pathway.getDisplayName());
            }
        }
        fu.close();
    }
    
    /**
     * Use this method to generate a map from FIs to Pathways.
     * @throws Exception
     */
    @Test
    public void generateFIsToPathwaysMap() throws Exception {
        Collection reactions = prepareReactions();
        Collection complexes = prepareComplexes();
        MySQLAdaptor dba = (MySQLAdaptor) getMySQLAdaptor();
        Collection<GKInstance> pathways = dba.fetchInstanceByAttribute(ReactomeJavaConstants.Pathway,
                                                              ReactomeJavaConstants.species,
                                                              "=",
                                                              48887);
        dba.loadInstanceAttributeValues(pathways, new String[]{
                ReactomeJavaConstants.hasEvent,
                ReactomeJavaConstants.dataSource
        });
        System.out.println("Total human pathways: " + pathways.size());
        // Remove disease pathway
        Long diseaseId = 1643685L;
        GKInstance disease = dba.fetchInstance(diseaseId);
        Set<GKInstance> diseasePathways = InstanceUtilities.getContainedEvents(disease);
        pathways.removeAll(diseasePathways);
        pathways.remove(disease);
        System.out.println("Remove disease pathways: " + pathways.size());
        // Get the directory
        String resultDir = FIConfiguration.getConfiguration().get("RESULT_DIR");
        String fileName = resultDir + File.separator + "ReactomeFIsToPathways_082216.txt";
        FileUtility fu = new FileUtility();
        fu.setOutput(fileName);
        // No header is required
        Set<GKInstance> interactors = new HashSet<GKInstance>();
        Set<String> interactions = new HashSet<String>();
        for (GKInstance pathway : pathways) {
            if(pathway.getAttributeValue(ReactomeJavaConstants.dataSource) != null)
                continue; // Use Reactome only
            System.out.println("Handling " + pathway);
            Set<GKInstance> pathwayEvents = InstanceUtilities.grepPathwayEventComponents(pathway);
            interactions.clear();
            for (GKInstance event : pathwayEvents) {
                if (!(event.getSchemClass().isa(ReactomeJavaConstants.ReactionlikeEvent)))
                    continue;
                interactors.clear();
                extractInteractorsFromReaction(event, interactors);
                generateInteractions(interactors, interactions, event, true);
            }
            Set<GKInstance> pathwayEntities = InstanceUtilities.grepPathwayParticipants(pathway);
            for (GKInstance entity : pathwayEntities) {
                if (!entity.getSchemClass().isa(ReactomeJavaConstants.Complex))
                    continue;
                interactors.clear();
                grepComplexComponents(entity, interactors);
                generateInteractions(interactors, interactions, entity, true);
            }
            if (interactions.size() == 0)
                continue;
            for (String interaction : interactions) {
                // Use this format so that we can use gene set based classes
                interaction = interaction.replace("\t", ",");
                fu.printLine(interaction + "\t" +  pathway.getDisplayName());
            }
        }
        fu.close();
    }
    
    public Set<String> extractInteractionSet() throws Exception {
        Collection reactions = prepareReactions();
        Collection complexes = prepareComplexes();
        return extractInteractionSet(reactions, 
                                     complexes);
    }

    private Set<String> extractInteractionSet(Collection reactions,
                                              Collection complexes) throws Exception {
        Set<String> interactions = new HashSet<String>();
        GKInstance rxn = null;
        Set<GKInstance> interactors = new HashSet<GKInstance>();
        long time1 = System.currentTimeMillis();
        for (Iterator it = reactions.iterator(); it.hasNext();) {
            rxn = (GKInstance) it.next();
            //System.out.println("Reaction: " + c++);
            extractInteractorsFromReaction(rxn, interactors);
            generateInteractions(interactors, interactions, rxn);
            interactors.clear();
        }
        System.out.println("Total interactions from reactions: " + interactions.size());
        if (!excludeComplex) {
            GKInstance complex = null;
            for (Iterator it = complexes.iterator(); it.hasNext();) {
                complex = (GKInstance) it.next();
                //System.out.println("Complex: " + c++ + " " + complex.getDBID());
                interactors.clear();
                grepComplexComponents(complex, interactors);
                // No need
                //if (interactors.size() > 10)
                //    continue; // cutoff set manually
                generateInteractions(interactors, interactions, complex);
            }
        }
        long time2 = System.currentTimeMillis();
        System.out.println("Time for looping: " + (time2 - time1));
        System.out.println("Total interactions from Reactome: " + interactions.size());
        return interactions;
    }
    
    @SuppressWarnings("rawtypes")
    public void grepComplexComponents(GKInstance complex, Set<GKInstance> interactors) throws Exception {
        Set<GKInstance> current = new HashSet<GKInstance>();
        current.add(complex);
        Set<GKInstance> next = new HashSet<GKInstance>();
        while (current.size() > 0) {
            for (Iterator it = current.iterator(); it.hasNext();) {
                GKInstance tmp = (GKInstance) it.next();
                List components = tmp.getAttributeValuesList(ReactomeJavaConstants.hasComponent);
                if (components == null || components.size() == 0)
                    continue;
                for (Iterator it1 = components.iterator(); it1.hasNext();) {
                    GKInstance tmp1 = (GKInstance) it1.next();
                    if (tmp1.getSchemClass().isa(ReactomeJavaConstants.EntityWithAccessionedSequence))
                        interactors.add(tmp1);
                    else if (tmp1.getSchemClass().isa(ReactomeJavaConstants.EntitySet))
                        interactors.add(tmp1);
                    else if (tmp1.getSchemClass().isa(ReactomeJavaConstants.Complex))
                        next.add(tmp1);
                }
            }
            current.clear();
            current.addAll(next);
            next.clear();
        }
    }
    
    protected void generateInteractions(Set<GKInstance> interactors, 
                                        Set<String> interactions,
                                        GKInstance source) throws Exception {
        generateInteractions(interactors, interactions, source, false);
    }
    
    private void generateInteractions(Set<GKInstance> interactors, 
                                      Set<String> interactions,
                                      GKInstance source,
                                      boolean useGeneNames) throws Exception {
        List<GKInstance> list = new ArrayList<GKInstance>(interactors);
        int size = list.size();
        for (int i = 0; i < size - 1; i++) {
            GKInstance interactor1 = list.get(i);
            for (int j = i + 1; j < size; j++) {
                GKInstance interactor2 = list.get(j);
                generateInteractions(interactor1, 
                                     interactor2,
                                     interactions,
                                     source,
                                     useGeneNames);
            }
        }
    }
    
    /**
     * Generate interactions from the passed interactions into interactions as strings.
     * @param interactors
     * @param interactions
     * @param source
     * @throws Exception
     */
    public void generateInteractionsWithDBNames(Set<GKInstance> interactors,
                                                Set<String> interactions,
                                                GKInstance source) throws Exception {
        List<GKInstance> list = new ArrayList<GKInstance>(interactors);
        int size = list.size();
        for (int i = 0; i < size - 1; i++) {
            GKInstance interactor1 = list.get(i);
            for (int j = i + 1; j < size; j++) {
                GKInstance interactor2 = list.get(j);
                generateInteractionsWithDBNames(interactor1, interactor2, interactions, source);
            }
        }        
    }
    
    private void generateInteractionsWithDBNames(GKInstance interactor1,
                                                 GKInstance interactor2,
                                                 Set<String> interactions,
                                                 GKInstance source) throws Exception {
        Set<GKInstance> refPepSeqs1 = grepRefPepSeqs(interactor1);
        if (refPepSeqs1.size() == 0)
            return;
        Set<GKInstance> refPepSeqs2 = grepRefPepSeqs(interactor2);
        if (refPepSeqs1.size() == 0)
            return;
        // Permutate members in these two sets
        int comp = 0;
        String pair = null;
        for (GKInstance ref1 : refPepSeqs1) {
            String uid1 = (String) ref1.getAttributeValue(ReactomeJavaConstants.identifier);
            if (uid1 == null)
                continue;
            String dbName1 = null;
            GKInstance db1 = (GKInstance) ref1.getAttributeValue(ReactomeJavaConstants.referenceDatabase);
            if (db1 != null)
                dbName1 = db1.getDisplayName();
            else
                dbName1 = "unknown";
            for (GKInstance ref2 : refPepSeqs2) {
                String uid2 = (String) ref2.getAttributeValue(ReactomeJavaConstants.identifier);
                if (uid2 == null)
                    continue;
                String dbName2 = null;
                GKInstance db2 = (GKInstance) ref2.getAttributeValue(ReactomeJavaConstants.referenceDatabase);
                if (db2 != null)
                    dbName2 = db2.getDisplayName();
                else
                    dbName2 = "unknown";
                comp = uid1.compareTo(uid2);
                if (comp < 0)
                    pair = dbName1 + ":" + uid1 + "\t" + dbName2 + ":" + uid2;
                else if (comp > 0)
                    pair = dbName2 + ":" + uid2 + "\t" + dbName1 + ":" + uid1;
                if (pair != null) {
                    interactions.add(pair); //exclude self interaction
                }
            }
        }
    }
    
    private void generateInteractions(GKInstance interactor1, 
                                      GKInstance interactor2,
                                      Set<String> interactions,
                                      GKInstance source,
                                      boolean useGeneNames) throws Exception {
        if (excludeComplex) {
            if (interactor1.getSchemClass().isa(ReactomeJavaConstants.Complex) ||
                interactor2.getSchemClass().isa(ReactomeJavaConstants.Complex))
                return;
        }
        Set<GKInstance> refPepSeqs1 = grepRefPepSeqs(interactor1);
        Set<GKInstance> refPepSeqs2 = grepRefPepSeqs(interactor2);
        generateFIs(refPepSeqs1, 
                    refPepSeqs2,
                    interactions,
                    useGeneNames);
    }

    protected void generateFIs(Set<GKInstance> refPepSeqs1, 
                               Set<GKInstance> refPepSeqs2, 
                               Set<String> interactions,
                               boolean useGeneNames) throws InvalidAttributeException, Exception {
        if (refPepSeqs1.size() == 0 || refPepSeqs2.size() == 0)
            return;
        // Permutate members in these two sets
        int comp = 0;
        String pair = null;
        for (GKInstance ref1 : refPepSeqs1) {
            String uid1 = null;
            if (useGeneNames && ref1.getSchemClass().isValidAttribute(ReactomeJavaConstants.geneName))
                uid1 = (String) ref1.getAttributeValue(ReactomeJavaConstants.geneName);
            else
                uid1 = (String) ref1.getAttributeValue(ReactomeJavaConstants.identifier);
            if (uid1 == null)
                continue;
            for (GKInstance ref2 : refPepSeqs2) {
                String uid2 = null;
                if (useGeneNames && ref2.getSchemClass().isValidAttribute(ReactomeJavaConstants.geneName))
                    uid2 = (String) ref2.getAttributeValue(ReactomeJavaConstants.geneName);
                else
                    uid2 = (String) ref2.getAttributeValue(ReactomeJavaConstants.identifier);
                if (uid2 == null)
                    continue;
                comp = uid1.compareTo(uid2);
                if (comp < 0)
                    pair = uid1 + "\t" + uid2;
                else if (comp > 0)
                    pair = uid2 + "\t" + uid1;
                if (pair != null) {
                    interactions.add(pair); //exclude self interaction
                    // Used for debugging
                    //if (pair.equals("O95405 P01270"))
                    //    System.out.println(pair + " < " + source.getDBID());
                }
            }
        }
    }
    
    protected Set<GKInstance> grepRefPepSeqs(GKInstance interactor) throws Exception {
        return topicHelper.grepRefPepSeqs(interactor);
    }
    
    protected Collection prepareComplexes() throws Exception {
        MySQLAdaptor dba = (MySQLAdaptor) getMySQLAdaptor();
        GKInstance homosapiens = dba.fetchInstance(48887L);
        Collection complexes = dba.fetchInstanceByAttribute(ReactomeJavaConstants.Complex,
                                                            ReactomeJavaConstants.species,
                                                            "=",
                                                            homosapiens);
        // Change made on January 31, 2017. Previously all human complexes are used,
        // which should be regarded as a bug. Do the following for filtering.
        for (Iterator it = complexes.iterator(); it.hasNext();) {
            GKInstance complex = (GKInstance) it.next();
            GKInstance dataSource = (GKInstance) complex.getAttributeValue(ReactomeJavaConstants.dataSource);
            if (dataSource != null)
                it.remove();
        }
        SchemaClass cls = dba.getSchema().getClassByName(ReactomeJavaConstants.Complex);
        SchemaAttribute att = cls.getAttribute(ReactomeJavaConstants.hasComponent);
        dba.loadInstanceAttributeValues(complexes, att);
        return complexes;
    }
    
    protected Collection prepareReactions() throws Exception {
        // Load necessary attributes
        MySQLAdaptor dba = (MySQLAdaptor) getMySQLAdaptor();
        // Load all reactions for analyzed
        GKInstance homosapiens = dba.fetchInstance(48887L);
        Collection reactions = null;
        SchemaClass reactionCls = null;
        // Adjust for new schema
        if (dba.getSchema().isValidClass(ReactomeJavaConstants.ReactionlikeEvent)) {
            reactions = dba.fetchInstanceByAttribute(ReactomeJavaConstants.ReactionlikeEvent,
                                                    ReactomeJavaConstants.species,
                                                    "=",
                                                    homosapiens);
            reactionCls = dba.getSchema().getClassByName(ReactomeJavaConstants.ReactionlikeEvent);
        }
        else {
            reactions = dba.fetchInstanceByAttribute(ReactomeJavaConstants.Reaction,
                                                     ReactomeJavaConstants.species,
                                                     "=",
                                                     homosapiens);
            reactionCls = dba.getSchema().getClassByName(ReactomeJavaConstants.Reaction);
        }
        // Need a little bit filtering for Reactome reactions only
        for (Iterator it = reactions.iterator(); it.hasNext();) {
            GKInstance rxt = (GKInstance) it.next();
            GKInstance dataSource = (GKInstance) rxt.getAttributeValue(ReactomeJavaConstants.dataSource);
            if (dataSource != null)
                it.remove();
        }
        Collection cas = dba.fetchInstancesByClass(ReactomeJavaConstants.CatalystActivity);
        Collection regulations = dba.fetchInstancesByClass(ReactomeJavaConstants.Regulation);
        Collection entities = dba.fetchInstanceByAttribute(ReactomeJavaConstants.EntityWithAccessionedSequence,
                                                           ReactomeJavaConstants.species,
                                                           "=",
                                                           homosapiens);
        SchemaAttribute att = reactionCls.getAttribute(ReactomeJavaConstants.input);
        dba.loadInstanceAttributeValues(reactions, att);
        att = reactionCls.getAttribute(ReactomeJavaConstants.output);
        dba.loadInstanceAttributeValues(reactions, att);
        att = reactionCls.getAttribute(ReactomeJavaConstants.catalystActivity);
        dba.loadInstanceAttributeValues(reactions, att);
        att = reactionCls.getAttribute(ReactomeJavaConstants.regulatedBy);
        dba.loadInstanceAttributeValues(reactions, att);
        reactionCls = dba.getSchema().getClassByName(ReactomeJavaConstants.CatalystActivity);
        att = reactionCls.getAttribute(ReactomeJavaConstants.physicalEntity);
        dba.loadInstanceAttributeValues(cas, att);
        reactionCls = dba.getSchema().getClassByName(ReactomeJavaConstants.Regulation);
//        att = reactionCls.getAttribute(ReactomeJavaConstants.regulatedEntity);
//        dba.loadInstanceAttributeValues(regulations, att);
        att = reactionCls.getAttribute(ReactomeJavaConstants.regulator);
        dba.loadInstanceAttributeValues(regulations, att);
        reactionCls = dba.getSchema().getClassByName(ReactomeJavaConstants.EntityWithAccessionedSequence);
        att = reactionCls.getAttribute(ReactomeJavaConstants.referenceEntity);
        dba.loadInstanceAttributeValues(entities, att);
        return reactions;
    }
    
    public void findReactionGraphComponents() throws Exception {
        Collection reactions = loadReactions();
        final Map<GKInstance, Set<GKInstance>> componentMap = new HashMap<GKInstance, Set<GKInstance>>();
        GKInstance rxn = null;
        Set<String> preTypes = new HashSet<String>();
        for (Iterator it = reactions.iterator(); it.hasNext();) {
            rxn = (GKInstance) it.next();
            Set<GKInstance> component = componentMap.get(rxn);
            if (component == null) {
                component = new HashSet<GKInstance>();
                component.add(rxn);
                componentMap.put(rxn, component);
            }
            List preceding = rxn.getAttributeValuesList(ReactomeJavaConstants.precedingEvent);
            if (preceding != null) {
                for (Iterator it1 = preceding.iterator(); it1.hasNext();) {
                    GKInstance preRxn = (GKInstance) it1.next();
                    preTypes.add(preRxn.getSchemClass().getName());
                    component.add(preRxn);
                    componentMap.put(preRxn, component);
                }
            }
        }
        Set<Set<GKInstance>> components = new HashSet<Set<GKInstance>>(componentMap.values());
        List<Set<GKInstance>> list = new ArrayList<Set<GKInstance>>(components);
        Collections.sort(list, new Comparator<Set<GKInstance>>() {
            public int compare(Set<GKInstance> set1, Set<GKInstance> set2) {
                return set2.size() - set1.size();
            }
        });
        System.out.printf("Total Components: %d -> %d%n", componentMap.size(), components.size());
        // Print out the first ten
        for (int i = 0; i < 10; i++) {
            Set<GKInstance> set = list.get(i);
            System.out.println(i + ": " + set.size());
        }
        System.out.println("Preceding Type: " + preTypes);
    }
    
    public void checkNoInteractionInReactome() throws IOException {
        String noIntFileName = "results/NoInteractionsForTrain.txt";
        String reactomeIntFileName = "results/ReactomeInteractions.txt";
        Set<String> noIntSet = new HashSet<String>();
        Set<String> rIntSet = new HashSet<String>();
        FileUtility fu = new FileUtility();
        fu.setInput(noIntFileName);
        String line = null;
        while ((line = fu.readLine()) != null) {
            noIntSet.add(line);
        }
        fu.close();
        fu.setInput(reactomeIntFileName);
        while ((line = fu.readLine()) != null) {
            rIntSet.add(line);
        }
        fu.close();
        System.out.println("Reactome Interactions: " + rIntSet.size());
        int size1 = noIntSet.size();
        System.out.println("No Interactions: " + size1);
        noIntSet.removeAll(rIntSet);
        int size2 = noIntSet.size();
        System.out.println("Reactome in No Interactions: " + (size1 - size2));
    }
    
    public void countTotalUsedUniProtIDs() throws Exception {
        // Fetch human only
        MySQLAdaptor dba = (MySQLAdaptor) getMySQLAdaptor();
        Collection ews = dba.fetchInstanceByAttribute(ReactomeJavaConstants.EntityWithAccessionedSequence,
                                                      ReactomeJavaConstants.species,
                                                      "=",
                                                      48887);
        SchemaClass cls = dba.getSchema().getClassByName(ReactomeJavaConstants.EntityWithAccessionedSequence);
        SchemaAttribute att = cls.getAttribute(ReactomeJavaConstants.referenceEntity);
        dba.loadInstanceAttributeValues(ews, att);
        Set<String> ids = new HashSet<String>();
        for (Iterator it = ews.iterator(); it.hasNext();) {
            GKInstance ew = (GKInstance) it.next();
            GKInstance ref = (GKInstance) ew.getAttributeValue(ReactomeJavaConstants.referenceEntity);
            if (ref == null)
                continue;
            String id = (String) ref.getAttributeValue(ReactomeJavaConstants.identifier);
            if (id != null)
                ids.add(id);
        }
        System.out.println("Total IDs used in gk_central: " + ids.size());
    }

    public GKInstance getDataSource() throws Exception {
        if (dataSourceId == null)
            return null;
        // GKInstance for dataSource pantherdb
        GKInstance dataSource = null;
        PersistenceAdaptor adaptor = getMySQLAdaptor();
        if (adaptor instanceof MySQLAdaptor) 
            dataSource = ((MySQLAdaptor)adaptor).fetchInstance(dataSourceId);
        else
            dataSource = ((XMLFileAdaptor)adaptor).fetchInstance(dataSourceId);
        return dataSource;
    }
    
    /**
     * @throws Exception
     */
    @Test
    public void countProteinsInPathways() throws Exception {
        // Top level pathways
        MySQLAdaptor dba = (MySQLAdaptor) getMySQLAdaptor ();
        Collection events = dba.fetchInstancesByClass(ReactomeJavaConstants.Event);
        List topLevelPathways = InstanceUtilities.grepTopLevelEvents(events);
        List<GKInstance> topics = getTopics();
        // Check proteins in each topics
        Set<String> totals = new HashSet<String>();
        for (GKInstance topic : topics) {
            Set<String> ids = grepIDsFromTopic(topic);
            System.out.println(topic.getDisplayName() + "(" + topic.getDBID() + "): " + ids.size());
            if (!topLevelPathways.contains(topic))
                System.out.println("\tNot a toplevel!");
            totals.addAll(ids);
            if (ids.size() > 200) {
                // check the subpathways
                List components = null;
                if (topic.getSchemClass().isValidAttribute(ReactomeJavaConstants.hasComponent)) 
                    components = topic.getAttributeValuesList(ReactomeJavaConstants.hasComponent);
                else if (topic.getSchemClass().isValidAttribute(ReactomeJavaConstants.hasEvent))
                    components = topic.getAttributeValuesList(ReactomeJavaConstants.hasEvent);
                else if (topic.getSchemClass().isValidAttribute(ReactomeJavaConstants.hasSpecialisedForm))
                    components = topic.getAttributeValuesList(ReactomeJavaConstants.hasSpecialisedForm);
                else if (topic.getSchemClass().isValidAttribute(ReactomeJavaConstants.hasMember))
                    components = topic.getAttributeValuesList(ReactomeJavaConstants.hasMember);
                if (components == null || components.size() == 0)
                    continue;
                for (Iterator it = components.iterator(); it.hasNext();) {
                    GKInstance comp = (GKInstance) it.next();
                    Set<String> subIds = grepIDsFromTopic(comp);
                    System.out.println("\t" + comp.getDisplayName() + "(" + comp.getDBID() + "): " + subIds.size());
                }
            }
        }
        System.out.println("Total Ids: " + totals.size());
        // Special handling on Signaling Pathways
        GKInstance signalingPathway = dba.fetchInstance(162582L);
        List components = signalingPathway.getAttributeValuesList(ReactomeJavaConstants.hasComponent);
        for (Iterator it = components.iterator(); it.hasNext();) {
            GKInstance comp = (GKInstance) it.next();
            Collection referrers = comp.getReferers(ReactomeJavaConstants.hasComponent);
            if (referrers.size() > 1) {
                System.out.println("\t" + comp.getDisplayName() + "(" + comp.getDBID() + ") contained by others");
                continue;
            }
            referrers = comp.getReferers(ReactomeJavaConstants.hasMember);
            if (referrers != null && referrers.size() > 0) {
                System.out.println("\t" + comp.getDisplayName() + "(" + comp.getDBID() + ") contained by others");
                continue;
            }
            referrers = comp.getReferers(ReactomeJavaConstants.hasSpecialisedForm);
            if (referrers != null && referrers.size() > 0) {
                System.out.println("\t" + comp.getDisplayName() + "(" + comp.getDBID() + ") contained by others");
                continue;
            }
        }
    }
    
    /**
     * Generate a list of pathways from the Reactome pathway hierarchy tree so that they can be used
     * in pathway enrichment analysis. The generation is done in two steps:
     * 1). Check item in the FrontPageItem list. If a pathway item has less than 200 proteins, that pathway
     * is listed in the pathway list.
     * 2). Otherwise, a front page pathway's sub-pathways are listed in the list regardless of their sizes if 
     * the size of subpathways are less than 300. If a sub-pathway has size > 300, its sub-pathways are listed.
     * 3). Note: The cutoff values 200 and 300 is chosen rather arbitrary. The reason why 300 is chosen for the
     * second level is based on assumption that the lower level pathways should be more like functional units
     * than upper level pathways.
     * @throws Exception
     * @deprecated: Use method in PathwayGeneSetGenerator instead.
     */
    @Test
    @Deprecated
    public void generateListOfPathways() throws Exception {
        Set<GKInstance> pathwaySet = new HashSet<GKInstance>();
        MySQLAdaptor releasedDBA = (MySQLAdaptor) getMySQLAdaptor();
        Collection frontPages = releasedDBA.fetchInstancesByClass(ReactomeJavaConstants.FrontPage);
        List<GKInstance> bigTopics = new ArrayList<GKInstance>();
        for (Iterator it = frontPages.iterator(); it.hasNext();) {
            GKInstance frontPage = (GKInstance) it.next();
            List items = frontPage.getAttributeValuesList(ReactomeJavaConstants.frontPageItem);
            for (Iterator it1 = items.iterator(); it1.hasNext();) {
                GKInstance topic = (GKInstance) it1.next();
                // Make sure this is a human pathway
                GKInstance species = (GKInstance) topic.getAttributeValue(ReactomeJavaConstants.species);
                if (!species.getDBID().equals(48887L))
                    continue;
                Set<String> ids = grepIDsFromTopic(topic);
                if (ids.size() > 200) // This is arbitrary
                    bigTopics.add(topic);
                else
                    pathwaySet.add(topic);
            }
        }
        Set<GKInstance> next = new HashSet<GKInstance>();
        for (int i = 0; i < 2; i++) { // Run two levels only
            // Have to split the big topics
            for (GKInstance topic : bigTopics) {
                List comps = null;
                if (topic.getSchemClass().isValidAttribute(ReactomeJavaConstants.hasEvent))
                    comps = topic.getAttributeValuesList(ReactomeJavaConstants.hasEvent);
                else if (topic.getSchemClass().isValidAttribute(ReactomeJavaConstants.hasSpecialisedForm))
                    comps = topic.getAttributeValuesList(ReactomeJavaConstants.hasSpecialisedForm);
                else if (topic.getSchemClass().isValidAttribute(ReactomeJavaConstants.hasMember))
                    comps = topic.getAttributeValuesList(ReactomeJavaConstants.hasMember);
                if (comps == null || comps.size() == 0) {
                    pathwaySet.add(topic);
                    continue;
                }
                // If there is any reaction in the pathway, it should not split any more
                boolean isAdded = false;
                for (Object obj : comps) {
                    GKInstance subEvent = (GKInstance) obj;
                    if (subEvent.getSchemClass().isa(ReactomeJavaConstants.ReactionlikeEvent)) {
                        pathwaySet.add(topic);
                        isAdded = true;
                        break;
                    }
                }
                if (isAdded)
                    continue;
                for (Iterator<?> it = comps.iterator(); it.hasNext();) {
                    GKInstance sub = (GKInstance) it.next();
                    if (i == 1)
                        pathwaySet.add(sub);
                    else {
                        //if (sub.getDBID().equals(163359L)) 
                        //    continue; // Escape Glucagon signaling in metabolic regulation(163359)
                        //              // This pathway has been included by Integration of pathways involved in energy metabolism(163685)
                        // Check sub-pathway size
                        Set<String> ids = grepIDsFromTopic(sub);
                        if (ids.size() > 300) 
                            next.add(sub);
                        else
                            pathwaySet.add(sub);
                    }
                }
            }
            bigTopics.clear();
            bigTopics.addAll(next);
        }
        if (true) {
            pathwaySet.forEach(pathway -> System.out.println(pathway.getDBID() + "\t" + pathway.getDisplayName()));
            return;
        }
        // Want to sort it before output
        List<GKInstance> list = new ArrayList<GKInstance>(pathwaySet);
        InstanceUtilities.sortInstances(list);
        String fileName = FIConfiguration.getConfiguration().get("REACTOME_PATHWAYS");
        FileUtility fu = new FileUtility();
        fu.setOutput(fileName);
        for (GKInstance pathway : list)
            fu.printLine(pathway.getDBID() + "\t" + pathway.getDisplayName());
        fu.close();
        System.out.println("Total Pathways: " + list.size());
    }
    
    public Map<String, String> getUniProtToGeneMap(MySQLAdaptor dba) throws Exception {
        Collection<GKInstance> refSeqs = dba.fetchInstanceByAttribute(ReactomeJavaConstants.ReferenceGeneProduct,
                                                                      ReactomeJavaConstants.species,
                                                                      "=",
                                                                      48887L);
        dba.loadInstanceAttributeValues(refSeqs, new String[]{
                ReactomeJavaConstants.identifier,
                ReactomeJavaConstants.geneName
        });
        Map<String, String> idToGene = new HashMap<String, String>();
        for (GKInstance inst : refSeqs) {
            String id = (String) inst.getAttributeValue(ReactomeJavaConstants.identifier);
            String gene = (String) inst.getAttributeValue(ReactomeJavaConstants.geneName);
            idToGene.put(id, gene);
        }
        return idToGene;
    }
    
    /**
     * This method is used to create a map from synonyms to gene symbol. Basically it checks attributes
     * in the geneName slot of ReferenceGeneProduct. The original gene names are quite messy. A primary gene
     * name may be used as synonym in other primary gene name. For example, PCM1 is a synonym of MBD1. Therefore,
     * the mapping is controlled to avoid a primary gene name is mapped to another primary gene name. See the code 
     * note below.
     * @param dba
     * @return
     * @throws Exception
     */
    public Map<String, String> getSynonymToGeneName(MySQLAdaptor dba) throws Exception {
        Collection<GKInstance> refSeqs = dba.fetchInstanceByAttribute(ReactomeJavaConstants.ReferenceGeneProduct,
                                                                      ReactomeJavaConstants.species,
                                                                      "=",
                                                                      48887L);
        dba.loadInstanceAttributeValues(refSeqs, new String[]{
                ReactomeJavaConstants.geneName
        });
        Map<String, String> synonymToGene = new HashMap<String, String>();
        // Get the primary gene names first
        Set<String> primaryGeneNames = new HashSet<>();
        for (GKInstance inst : refSeqs) {
            if (inst.getSchemClass().isa(ReactomeJavaConstants.ReferenceIsoform))
                continue; // Just check ReferenceGeneProduct only
            List<String> geneNames = inst.getAttributeValuesList(ReactomeJavaConstants.geneName);
            if (geneNames == null || geneNames.size() == 0)
                continue;
            String gene = geneNames.get(0); // The first should be the current version of symbol
            primaryGeneNames.add(gene);
        }
        // To make it consistent, we will sort instances
        List<GKInstance> refSeqsList = refSeqs.stream()
                .filter(i -> !i.getSchemClass().isa(ReactomeJavaConstants.ReferenceIsoform))
                .sorted((i1, i2) -> i1.getDBID().compareTo(i2.getDBID()))
                .collect(Collectors.toList());
        // Now do the synonym mapping
        for (GKInstance inst : refSeqsList) {
            if (inst.getSchemClass().isa(ReactomeJavaConstants.ReferenceIsoform))
                continue; // Just check ReferenceGeneProduct only
            List<String> geneNames = inst.getAttributeValuesList(ReactomeJavaConstants.geneName);
            if (geneNames == null || geneNames.size() == 0)
                continue;
            String gene = geneNames.get(0); // The first should be the current version of symbol
            for (String name : geneNames) {
                // Don't map the primary gene name for the time being so that if it is used as a synonym
                // it will not be pushed into the map
                if (primaryGeneNames.contains(name))
                    continue;
                // Because one gene may be mapped to more than on UniProt, 
                // a synonym or a gene may be in the map already. Since we have sorted these instances,
                // therefore, the first primary gene will be used in the ordered list. This may not be correct
                // for some genes. But it should make a consistent result.
                if (synonymToGene.containsKey(name)) {
//                    logger.warn(name + " has been mapped to a primary gene name earlier.");
                    continue;
                }
                synonymToGene.put(name, gene);
            }
        }
        // Now push all primary gene names by themselves
        primaryGeneNames.forEach(name -> synonymToGene.put(name, name));
        return synonymToGene;
    }
}
