/*
 * Created on Apr 25, 2007
 *
 */
package org.reactome.fi;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.gk.model.GKInstance;
import org.junit.Test;
import org.reactome.data.ProteinIdFilters;
import org.reactome.data.ReactomeAnalyzer;
import org.reactome.data.UniProtAnalyzer;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.FeatureChecker;
import org.reactome.fi.util.FileUtility;
import org.reactome.fi.util.InteractionUtilities;
import org.reactome.fi.util.PositiveChecker;
import org.reactome.fi.util.Value;
import org.reactome.hibernate.HibernateFIReader;
import org.reactome.tred.TREDAnalyzer;
import org.reactome.weka.FeatureHandlerForV3;

/**
 * This class is used to handle the interaction file generated from the hibernate adaptor.
 * @author guanming
 *
 */
public class FIFileAnalyzer {
    private final String pathwayFIFile = FIConfiguration.getConfiguration().get("RESULT_DIR") + "PathwayFIs040909.txt";
    private final String predictedFIFile = FIConfiguration.getConfiguration().get("PREDICTED_FI_FILE");
    
    private FileUtility fu;
    
    public FIFileAnalyzer() {
        fu = new FileUtility();
    }
    
    @Test
    public void checkKEGGFIsInOthers() throws Exception {
        Set<String> predictedFIs = fu.loadInteractions(FIConfiguration.getConfiguration().get("GENE_FI_PREDICTED_FILE_NAME"));
        System.out.println("Total predicted FIs in genes: " + predictedFIs.size());
        Set<String> pathwayFIs = fu.loadInteractions(FIConfiguration.getConfiguration().get("GENE_FI_PATHWAY_FILE_NAME"));
        System.out.println("Total pathway FIs in genes: " + pathwayFIs.size());
        Set<String> keggFIs = fu.loadInteractions(FIConfiguration.getConfiguration().get("RESULT_DIR") + "/../2022/FIs_KEGG Pathway.txt");
        System.out.println("Total KEGG FIs in UniProt: " + keggFIs.size());
        // Do a mapping
        HibernateFIReader fiReader = new HibernateFIReader();
        Map<String, String> id2name = fiReader.generateAccessionToProteinNames();
        Set<String> keggFIsInGenes = new HashSet<>();
        for (String fi : keggFIs) {
            String[] tokens = fi.split("\t");
            String gene1 = id2name.get(tokens[0]);
            String gene2 = id2name.get(tokens[1]);
            if (gene1 == null || gene2 == null || gene1.equals(gene2))
                continue;
            String geneFI = InteractionUtilities.generateFIFromGene(gene1, gene2);
            if (geneFI == null)
                continue;
            keggFIsInGenes.add(geneFI);
        }
        System.out.println("Kegg FIs in genes: " + keggFIsInGenes.size());
        Set<String> copy = new HashSet<>(keggFIsInGenes);
        copy.removeAll(pathwayFIs);
        System.out.println("After removing pathway fis: " + copy.size());
        copy.removeAll(predictedFIs);
        System.out.println("After removing predicted FIs: " + copy.size());
    }
    
    @Test
    public void checkUsedReactomePathwaysAndGenes() throws Exception {
        String fileName = "/Users/wug/git/FIVizWS_corews/src/main/webapp/WEB-INF/ProteinNameToReactomePathways_Rel_71_091720.txt";
        Set<String> genes = new HashSet<>();
        Set<String> pathways = new HashSet<>();
        Files.lines(Paths.get(fileName))
             .map(line -> line.split("\t"))
             .forEach(tokens -> {
                 genes.add(tokens[0]);
                 pathways.add(tokens[1]);
             });
        System.out.println("Total genes: " + genes.size());
        System.out.println("Total pathways: " + pathways.size());
    }
    
    @Test
    public void addHasPPIEvidenceFeature() throws Exception {
        String dirName = "/Users/gwu/git/Ogmios/results/";
        String inFile = dirName + "ProteinFIsInReactions_032017.txt";
        
        Set<String> fis = new HashSet<String>();
        fu.setInput(inFile);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            fis.add(tokens[0] + "\t" + tokens[1]);
        }
        fu.close();
        // Attach PPI feature
        FeatureHandlerForV3 featureHandler = new FeatureHandlerForV3();
        List<String> featureList = featureHandler.getFeatureList();
        Map<String, Value> fiToValue = featureHandler.convertPairsToValues(fis, true);
        
        String outFile = dirName + "ProteinFIsInReactionsWithPPIEvidence_032017.txt";
        fu.setInput(inFile);
        fu.setOutput(outFile);
        line = fu.readLine();
        fu.printLine(line + "\thasPPIEvidence");
        int totalPPITrue = 0;
        int totalFIs = 0;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            Value value = fiToValue.get(tokens[0] + "\t" + tokens[1]);
            Boolean posFeature = value.humanInteraction |
                    value.mousePPI |
                    value.dmePPI |
                    value.celPPI |
                    value.scePPI |
                    value.pfamDomainInt;
            fu.printLine(line + "\t" + posFeature);
            totalFIs ++;
            if (posFeature)
                totalPPITrue ++;
        }
        fu.close();
        System.out.println("Total FIs: " + totalFIs);
        System.out.println("Total FIs having PPI evidence: " + totalPPITrue);
    }
    
//    /**
//     * This method is used to remove the ZNF clique.
//     * @throws IOException
//     */
//    @Test
//    public void generateFIFileWithZNFRemoved() throws IOException {
//        Set<String> fis = fu.loadInteractions(FIConfiguration.getConfiguration().get("GENE_FI_FILE_NAME"));
//        Collection<String> znfClique = new GraphAnalyzer().searchZNFClique(FIConfiguration.getConfiguration().get("GENE_FI_FILE_NAME"));
////        System.out.println("Genes in the ZNF clique: " + znfClique.size());
////        System.out.println(znfClique);
//        Set<String> znfFIs = InteractionUtilities.getFIs(znfClique, fis);
//        fis.removeAll(znfFIs);
//        String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR + "FIsInGene_No_ZNF_042810.txt";
//        fu.saveInteractions(fis, fileName);
//    }
    
    /**
     * This method is used to check the connections to UBC, which has over 3,000 FIs
     * in the 2015 version of the Reactome FI network.
     * @throws IOException
     */
    @Test
    public void checkUBCConnections() throws IOException {
//        String fiFileName = FIConfiguration.getConfiguration().get("GENE_FI_FILE_NAME");
        String fiFileName = FIConfiguration.getConfiguration().get("GENE_FI_BIG_COMP_FILE_NAME");
        System.out.println("FI file name: " + fiFileName);
        Set<String> fis = fu.loadInteractions(fiFileName);
        Map<String, Set<String>> proteinToPartners = InteractionUtilities.generateProteinToPartners(fis);
        Set<String> ubcPartners = proteinToPartners.get("UBC");
        System.out.println("UBC: " + ubcPartners.size());
        
        // Generate a file without UBC
        int lastIndex = fiFileName.lastIndexOf(".");
        String fiFileName1 = fiFileName.substring(0, lastIndex) + "_NoUBC" + fiFileName.substring(lastIndex);
        fu.setOutput(fiFileName1);
        for (String fi : fis) {
            String[] tokens = fi.split("\t");
            if (tokens[0].equals("UBC") || tokens[1].equals("UBC"))
                continue;
            fu.printLine(fi);
        }
        fu.close();
    }
    
    /**
     * This method is used to check ZNF genes in the FI files.
     * @throws IOException
     */
    @Test
    public void checkZNFInFIs() throws IOException {
        Set<String> fis = fu.loadInteractions(FIConfiguration.getConfiguration().get("GENE_FI_FILE_NAME"));
        Set<String> genes = InteractionUtilities.grepIDsFromInteractions(fis);
        Set<String> znfs = new HashSet<String>();
        for (String gene : genes) {
            if (gene.startsWith("ZNF"))
                znfs.add(gene);
        }
        System.out.println("Total Genes: " + genes.size());
        System.out.println("Total ZNFs: " + znfs.size() + " (" + (double)znfs.size() / genes.size() + ")");
        Set<String> touchedZNFFIs = new HashSet<String>();
        Set<String> allZNFFIs = new HashSet<String>();
        for (String fi : fis) {
            int index = fi.indexOf("\t");
            String gene1 = fi.substring(0, index);
            String gene2 = fi.substring(index + 1);
            if (gene1.startsWith("ZNF") || gene2.startsWith("ZNF"))
                touchedZNFFIs.add(fi);
            if (gene1.startsWith("ZNF") && gene2.startsWith("ZNF"))
                allZNFFIs.add(fi);
        }
        System.out.println("Total FIs: " + fis.size());
        System.out.println("Total FIs having at least one ZNF: " + touchedZNFFIs.size() + " (" + (double)touchedZNFFIs.size() / fis.size() + ")");
        System.out.println("Total FIs having both ZNFs: " + allZNFFIs.size() + " (" + (double)allZNFFIs.size() / fis.size() + ")");
    }
    
    /**
     * Use this method to save a set of FIs into an order list.
     * @param fis
     * @param fileName
     * @throws IOException
     */
    public void saveFIInOrder(Set<String> fis,
                              String fileName) throws IOException {
        List<String> list = new ArrayList<String>(fis);
        Collections.sort(list);
        fu.setOutput(fileName);
        for (String fi : list)
            fu.printLine(fi);
        fu.close();
    }
    
    /**
     * Load all FIs in UniProt accession numbers.
     * @return
     * @throws IOException
     */
//    public Set<String> loadFIs() throws IOException {
//        Set<String> fis = new HashSet<String>();
//        fis.addAll(loadPathwayFIs());
//        fis.addAll(loadTFTargetInteractions());
//        fis.addAll(loadPredictedFIs());
//        return fis;
//        //return fu.loadInteractions(FIConfiguration.getConfiguration().get("INTERACTION_FILE_NAME);
//    }
    
    private Set<String> loadArtificalFIs() throws IOException {
        Set<String> fis = new HashSet<String>();
        fis.add("A B");
        fis.add("A C");
        fis.add("A D");
        fis.add("A G");
        fis.add("B D");
        fis.add("C D");
        fis.add("C F");
        fis.add("D E");
        fis.add("D G");
        fis.add("D H");
        fis.add("E F");
        fis.add("E H");
        fis.add("G H");
        return fis;
    }
    
    /**
     * This class is used to check the SwissProt coverage in the functional
     * interaction file.
     * @throws IOException
     */
    @Test
    public void checkSwissProtCoverage() throws Exception {
        UniProtAnalyzer uniProtAnalyzer = new UniProtAnalyzer();
        // Need to use this map in case some UniProt accession numbers
        // used in pathways are not the first one!
        Map<String, String> map = uniProtAnalyzer.loadSwissProtIDsMap();
        Set<String> mapIds = new HashSet<String>(map.values());
        // There are three files that need to be check
        FIConfiguration config = FIConfiguration.getConfiguration();
//        String[] fileNames = new String[] {
////                "FIInteractions73_021108_Pathway.txt",
////                "FIInteractions73_021108_PPI.txt",
////                "FIInteractions73_021108.txt"
////                "FI73_041408.txt"
//                config.get("")
//        };
//        for (String name : fileNames) {
//            checkSwissProtCoverage(map, 
//                                   mapIds, 
//                                   name);
//        }
        List<ReactomeAnalyzer> analyzers = ReactomeAnalyzer.getPathwayDBAnalyzers();
        Set<String> totalFIs = new HashSet<String>();
        for (ReactomeAnalyzer analyzer : analyzers) {
            String sourceName = getDataSourceName(analyzer);
            String fileName = config.get("RESULT_DIR") + "/FIs_" + sourceName + ".txt";
            checkSwissProtCoverage(map, mapIds, fileName);
//            if (sourceName.contains("ENCODE"))
//                continue;
            totalFIs.addAll(fu.loadInteractions(fileName));
        }
        System.out.println("Total Pathway FIs: ");
        checkSwissProtCoverage(map, mapIds, totalFIs);
        String fileName = config.get("PREDICTED_FI_FILE");
        Set<String> predictedFIs = fu.loadInteractions(fileName);
        System.out.println("\nTotal predicted FIs:");
        checkSwissProtCoverage(map, mapIds, predictedFIs);
        System.out.println("\nTotal FIs:");
        totalFIs.addAll(predictedFIs);
        checkSwissProtCoverage(map, mapIds, totalFIs);
    }

    private void checkSwissProtCoverage(Map<String, String> map,
                                        Set<String> mapIds, 
                                        String fileName) throws IOException {
        System.out.println("File: " + fileName);
        Set<String> interactions = fu.loadInteractions(fileName);
        checkSwissProtCoverage(map, mapIds, interactions);
    }

    private void checkSwissProtCoverage(Map<String, String> map,
                                        Set<String> mapIds,
                                        Set<String> interactions) {
        System.out.println("Total interactions: " + interactions.size());
        Set<String> totalIds = InteractionUtilities.grepIDsFromInteractions(interactions);
        System.out.println("Total IDs: " + totalIds.size());
        totalIds = removeSpliceIsoform(totalIds);
        System.out.println("Remove isoforms: " + totalIds.size());
        // 25205 is the total identifiers in HPRD and used as the total
        // gene numbers
        System.out.println("Total coverage: " + totalIds.size() / 25205.0);
        // Check in other way
        Set<String> tmpIds = new HashSet<String>();
        for (String id : totalIds) {
            String tmp = map.get(id);
            if (tmp != null)
                tmpIds.add(tmp);
        }
        System.out.println("Swiss Prot Ids in FIs: " + tmpIds.size());
        System.out.println("Swiss Coverage: " + tmpIds.size() + "/" + 
                           mapIds.size() + "=" + (double)tmpIds.size() / mapIds.size());
        System.out.println();
    }
    
    private Set<String> removeSpliceIsoform(Set<String> ids) {
        Set<String> rtn = new HashSet<String>();
        int index = 0;
        for (String id : ids) {
            index = id.indexOf("-");
            if (index > 0)
                rtn.add(id.substring(0, index));
            else
                rtn.add(id);
        }
        return rtn;
    }
    
//    public Set<String> loadInteractionIds() throws IOException {
//        Set<String> fis = loadFIs();
//        return InteractionUtilities.grepIDsFromInteractions(fis);
//    }
    
//    @Test
//    public void checkHumanIDsInFIs() throws IOException {
//        UniProtAnalyzer uniAnalyzer = new UniProtAnalyzer();
//        Map<String, String> uniProtIds = uniAnalyzer.loadUniProtIDsMap();
//        Set<String> uniSet = uniProtIds.keySet();
//        Set<String> fis = loadFIs();
//        int index = 0;
//        int totalFIs = 0;
//        int removeFIs = 0;
//        Set<String> idsFromHPRDOrNCBI = new HashSet<String>();
//        Set<String> removedIds = new HashSet<String>();
//        // Need to get alternative form out
//        for (String pair : fis) {
//            totalFIs ++;
//            index = pair.indexOf(" ");
//            String id1 = pair.substring(0, index);
//            String id2 = pair.substring(index + 1);
//            if (id1.matches("^[0-9]+")) {
//                idsFromHPRDOrNCBI.add(id1);
//            }
//            if (id2.matches("^[0-9]+"))
//                idsFromHPRDOrNCBI.add(id2);
//            index = id1.indexOf("-");
//            if (index > 0)
//                id1 = id1.substring(0, index);
//            index = id2.indexOf("-");
//            if (index > 0)
//                id2 = id2.substring(0, index);
//            // Check if id1 or id2 are numbers only: for NCBI or HPRD
//            if (!id1.matches("^[0-9]+") && !uniSet.contains(id1)) {
//                removedIds.add(id1);
//                removeFIs ++;
//                continue;
//            }
//            if (!id2.matches("^[0-9]+") && !uniSet.contains(id2)) {
//                removedIds.add(id2);
//                removeFIs ++;
//                continue;
//            }
//            // Otherwise, have to make sure they are from human UniProt IDs
//        }
//        System.out.println("Total FIs: " + totalFIs);
//        System.out.println("Remove FIs: " + removeFIs);
//        System.out.println("Total IDs from NCBI or HPRD: " + idsFromHPRDOrNCBI.size() + ": " + idsFromHPRDOrNCBI);
//        System.out.println("Removed Ids: " + removedIds.size());
//        for (String id : removedIds)
//            System.out.println(id);
//    }
    
//    public Map<String, Set<String>> loadIdToPartners() throws IOException {
//        Map<String, Set<String>> map = new HashMap<String, Set<String>>();
//        Set<String> fis = loadFIs();
//        return new BreadthFirstSearch().generateIdToPartnersMap(fis);
//    }
    
//    
//    @Test
//    public void analyzeNonPathwayIds() throws IOException {
//        Set<String> totalIds = loadInteractionIds();
//        TopicAnalyzer topicAnalyzer = new TopicAnalyzer();
//        Set<String> topicIds = topicAnalyzer.getTopicIds();
//        System.out.println("Total Ids: " + totalIds.size());
//        System.out.println("Pathway Ids: " + topicIds.size());
//        totalIds.removeAll(topicIds);
//        System.out.println("Non Pathway Ids: " + totalIds.size());
//        Map<String, Set<String>> idToPartners = loadIdToPartners();
//        int hop = 1;
//        while (totalIds.size() > 0) {
//            System.out.println("Hop: " + hop);
//            Set<String> hopAnnotated = checkNonPathwayIds(totalIds, 
//                                                          topicIds, 
//                                                          idToPartners);
//            if (hopAnnotated.size() == 0) {
//                // Cannot be annotated by hop anymore
//                System.out.println("Cannot annotated by hop!");
//                break;
//            }
//            System.out.println("Annotated by hopping: " + hopAnnotated.size());
//            totalIds.removeAll(hopAnnotated);
//            System.out.println("Not Annotated: " + totalIds.size());
//            topicIds.addAll(hopAnnotated);
//            hop ++;
//        }
//    }
    
    private Set<String> checkNonPathwayIds(Set<String> nonAnnotatedIds, 
                                           Set<String> topicIds, 
                                           Map<String, Set<String>> idToPartners) {
        Set<String> annotatedAfterHop = new HashSet<String>();
        for (Iterator<String> it = nonAnnotatedIds.iterator(); it.hasNext();) {
            String id = it.next();
            Set<String> partners = idToPartners.get(id);
            // Check if any partners is annotated
            for (String partner : partners) {
                if (topicIds.contains(partner)) {
                    annotatedAfterHop.add(id);
                    break;
                }
            }
        }
        return annotatedAfterHop;
    }
    
    /**
     * Generate a FI file in upper case
     * @throws IOException
     */
    @Test
    public void upperCaseForFIFile() throws IOException {
        String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR") + "FI73InGeneUpperCase_111208.txt";
        String inFile = FIConfiguration.getConfiguration().get("RESULT_DIR") + "FI73InGene_102908.txt";
        Set<String> fis = fu.loadInteractions(inFile);
        Set<String> newFis = new HashSet<String>();
        for (String fi : fis) {
            int index = fi.indexOf("\t");
            String name1 = fi.substring(0, index).toUpperCase();
            String name2 = fi.substring(index + 1).toUpperCase();
            int compare = name1.compareTo(name2);
            if (compare < 0)
                newFis.add(name1 + "\t" + name2);
            else if (compare > 0)
                newFis.add(name2 + "\t" + name1);
        }
        saveFIInOrder(newFis, fileName);
    }
    
    /**
     * This method is used to check the name case (upper or lower) usage for protein or
     * gene names in a FI file.
     * @throws IOException
     */
    @Test
    public void analyzeCaseInFIInNames() throws IOException {
        String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR") + "FI73InGene_102908.txt";
        Set<String> fis = fu.loadInteractions(fileName);
        Set<String> names = InteractionUtilities.grepIDsFromInteractions(fis);
        Map<String, Set<String>> upperToNames = new HashMap<String, Set<String>>();
        for (String name : names) {
            String upper = name.toUpperCase();
            Set<String> set = upperToNames.get(upper);
            if (set == null) {
                set = new HashSet<String>();
                upperToNames.put(upper, set);
            }
            set.add(name);
        }
        // Check more than one case
        System.out.println("Names having more than one cases:");
        for (String upper : upperToNames.keySet()) {
            Set<String> set = upperToNames.get(upper);
            if (set.size() > 1) {
                System.out.println(upper + ": " + set);
            }
        }
        System.out.println("\n\n");
        System.out.println("Names using lower cases:");
        for (String upper : upperToNames.keySet()) {
            Set<String> set = upperToNames.get(upper);
            if (set.size() == 1) {
                String name = set.iterator().next();
                if (!name.equals(upper))
                    System.out.println(upper + ": " + name);
            }
        }
        Map<String, Set<String>> proteinToPartners = InteractionUtilities.generateProteinToPartners(fis);
        // Check two cases have the same interactions
        System.out.println("\n\nCheck partners:");
        for (String upper : upperToNames.keySet()) {
            Set<String> set = upperToNames.get(upper);
            if (set.size() > 1) {
                System.out.println(upper);
                Iterator<String> it = set.iterator();
                String name1 = it.next();
                String name2 = it.next();
                Set<String> partners1 = proteinToPartners.get(name1);
                Set<String> partners2 = proteinToPartners.get(name2);
                if (partners1.equals(partners2))
                    System.out.println(name1 + ", " + name2 + " have the same partners!");
                else {
                    Set<String> shared = new HashSet<String>(partners1);
                    shared.retainAll(partners2);
                    partners1.removeAll(shared);
                    partners2.removeAll(shared);
                    System.out.println(name1 + ": " + partners1);
                    System.out.println(name2 + ": " + partners2);
                }
            }
        }
    }
    
    /**
     * This method is used to create files for pathway FIs.
     * @throws Exception
     */
    @Test
    public void dumpPathwayFIs() throws Exception {
        List<ReactomeAnalyzer> analyzers = ReactomeAnalyzer.getPathwayDBAnalyzers();
        ProteinIdFilters filters = new ProteinIdFilters();
        for (ReactomeAnalyzer analyzer : analyzers) {
            Set<String> fis = analyzer.extractInteractionSet();
            Set<String> normalized = filters.normalizeProteinPairs(fis);
            String sourceName = getDataSourceName(analyzer);
            System.out.println("Done data source: " + sourceName);
            String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR") + "/FIs_" + sourceName + ".txt";
            fu.saveInteractions(normalized, fileName);
            System.out.println();
        }
    }

    private String getDataSourceName(ReactomeAnalyzer analyzer)
            throws Exception {
        GKInstance dataSource = analyzer.getDataSource();
        String sourceName = null;
        if (dataSource == null)
            sourceName = "Reactome";
        else
            sourceName = dataSource.getDisplayName();
        return sourceName;
    }
    
    
    
    /**
     * Merge all pre-dumped pathway FI files into one: PathwayFIs????.txt.
     * @throws Exception
     */
    @Test
    public void generateOnePathwayFIFile() throws Exception {
        List<ReactomeAnalyzer> analyzers = ReactomeAnalyzer.getPathwayDBAnalyzers();
        Set<String> pathwayFIs = new HashSet<String>();
        for (ReactomeAnalyzer analyzer : analyzers) {
            String sourceName = getDataSourceName(analyzer);
            String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR") + "FIs_" + sourceName + ".txt";
            Set<String> fis = fu.loadInteractions(fileName);
            pathwayFIs.addAll(fis);
        }
        fu.saveInteractions(pathwayFIs, pathwayFIFile);
    }
    
    @Test
    public void checkTotalPathwayFIs() throws Exception {
        Set<String> allFIs = loadPathwayFIsFromFiles();
        System.out.println("Total pathway FIs: " + allFIs.size());
        Set<String> tredFIs = loadTFTargetInteractions();
        System.out.println("TRED FIs: " + tredFIs.size());
        allFIs.addAll(tredFIs);
        System.out.println("After merging: " + allFIs.size());
        ProteinAndInteractionCount counter = new ProteinAndInteractionCount();
        counter.countVsSwissProt(InteractionUtilities.grepIDsFromInteractions(allFIs));
    }
    
    @Test
    public void checkTotalFIsAndNamesInFile() throws IOException {
        String fileName = FIConfiguration.getConfiguration().get("GENE_FI_ANNOTATION_FILE_NAME");
        fu.setInput(fileName);
        Set<String> totalNames = new HashSet<String>();
        String line = fu.readLine();
        Set<String> fis = new HashSet<String>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            totalNames.add(tokens[0]);
            totalNames.add(tokens[1]);
            fis.add(tokens[0] + "\t" + tokens[1]);
        }
        fu.close();
        System.out.println("Total names: " + totalNames.size());
        System.out.println("Total FIs: " + fis.size());
    }
    
    /**
     * Load pathway FIs in UniProt ids.
     * @return
     * @throws IOException
     */
    @Deprecated
    public Set<String> loadPathwayFIs() throws Exception {
        return loadPathwayFIsFromFiles();
//        return fu.loadInteractions(pathwayFIFile);
    }
    
//    @Test
//    public void generateFIFiles() throws Exception {
//        Set<String> fis = loadPathwayAndTFTargetFIs();
//        // This file contains both pathways and TF/Target FIs.
//        fu.saveInteractions(fis, FIConfiguration.getConfiguration().get("RESULT_DIR") + "FIs_Pathway_043009.txt");
//        Set<String> predictedFIs = loadPredictedFIs();
//        fu.saveInteractions(predictedFIs, FIConfiguration.getConfiguration().get("RESULT_DIR") + "FIs_Predicted_043009.txt");
//        fis.addAll(predictedFIs);
//        fu.saveInteractions(fis, FIConfiguration.getConfiguration().get("RESULT_DIR") + "FIs_043009.txt");
//    }

    public Set<String> loadPathwayAndTFTargetFIs() throws Exception {
        // Pathway FIs
        Set<String> fis = loadPathwayFIsFromFiles();
        Set<String> tredFIs = loadTFTargetInteractions();
        fis.addAll(tredFIs);
        return fis;
    }
    
    /**
     * Load predicted FIs in UniProt ids.
     * @return
     * @throws IOException
     */
    public Set<String> loadPredictedFIs() throws IOException {
        return fu.loadInteractions(predictedFIFile);
    }
    
    /**
     * Load TF/Target interactions.
     * @return
     * @throws IOException
     */
    public Set<String> loadTFTargetInteractions() throws IOException {
        return fu.loadInteractions(FIConfiguration.getConfiguration().get("TRED_FI_FILE"));
        //        return fu.loadInteractions(FIConfiguration.getConfiguration().get("RESULT_DIR + "TREDInteractionsInUniProt.txt");
    }
    
    /**
     * This file is used to check the final FI network size.
     * @throws Exception
     */
    @Test
    public void checkFinalNetworkSize() throws Exception {
        ProteinAndInteractionCount counter = new ProteinAndInteractionCount();
        Set<String> pathwayFIs = loadPathwayFIs();
        System.out.println("Pathway FIs:");
        countFinalNetwork(counter, pathwayFIs);
        Set<String> predictedFIs = loadPredictedFIs();
        System.out.println("Predicted FIs:");
        countFinalNetwork(counter, predictedFIs);
        Set<String> tfTargetInteractions = new TREDAnalyzer().loadTFTargetInteractions();
        System.out.println("TF/Target interactions:");
        countFinalNetwork(counter, tfTargetInteractions);
        // Pathways and predicted
        predictedFIs.addAll(pathwayFIs);
        System.out.println("Pathway + Predicted FIs:");
        countFinalNetwork(counter, predictedFIs);
        // Pathways, predicted and TF/Targets
        predictedFIs.addAll(tfTargetInteractions);
        System.out.println("Pathway + Predicted + TF/Target FIs:");
        countFinalNetwork(counter, predictedFIs);
        // Count pathways and TF/Targets
        tfTargetInteractions.addAll(pathwayFIs);
        System.out.println("Pathway + TF/Target FIs:");
        countFinalNetwork(counter, tfTargetInteractions);
    }

    private void countFinalNetwork(ProteinAndInteractionCount counter,
                                   Set<String> pathwayFIs) throws IOException {
        Set<String> pathwayIds = InteractionUtilities.grepIDsFromInteractions(pathwayFIs);
        System.out.println("Total FIs: " + pathwayFIs.size() + " (" + pathwayIds.size() + ")");
        counter.countVsSwissProt(pathwayIds);
        System.out.println();
    }

    /**
     * Load a set of pre-generated normalized pathway FIs from files.
     * As of September, 2012, TF/Target interactions have been included
     * in pathway FIs.
     * @return
     * @throws IOException
     */
    public Set<String> loadPathwayFIsFromFiles() throws Exception {
        List<ReactomeAnalyzer> analyzers = ReactomeAnalyzer.getPathwayDBAnalyzers();
        FIConfiguration config = FIConfiguration.getConfiguration();
        Set<String> allFIs = new HashSet<String>();
        for (ReactomeAnalyzer analyzer : analyzers) {
            String sourceName = getDataSourceName(analyzer);
            String fileName = config.get("RESULT_DIR") + "/FIs_" + sourceName + ".txt";
            Set<String> fis = fu.loadInteractions(fileName);
            allFIs.addAll(fis);
        }
        return allFIs;
    }
    
    @Test
    public void testLoadPathwayFIsFromFiles() throws Exception {
        Set<String> pathwayFIs = loadPathwayFIsFromFiles();
        System.out.println("Total pathway FIs: " + pathwayFIs.size());
    }
    
    /**
     * This method is used to compare two versions of FIs to see if predicted FIs have
     * been confirmed by annotated FIs.
     * @throws Exception
     */
    @Test
    public void compareTwoVersionsOfFIs() throws Exception {
        // Load current version of FIs
        final Set<String> currentPathwayFIs = loadPathwayFIsFromFiles();
        System.out.println("Current pathway FIs: " + currentPathwayFIs.size());
        Set<String> currentPredictedFIs = loadPredictedFIs();
        System.out.println("Current predicted FIs: " + currentPredictedFIs.size());
        Set<String> currentFIs = new HashSet<String>(currentPathwayFIs);
        currentFIs.addAll(currentPredictedFIs);
        System.out.println("Total: " + currentFIs.size());
        // Load previous version of FIs
        String dirName = "/Users/gwu/Documents/EclipseWorkspace/caBigR3/results/v3/";
        Set<String> previousPathwayFIs = fu.loadInteractions(dirName + "FIs_Pathway_043009.txt");
        previousPathwayFIs = replaceSpaceWithTab(previousPathwayFIs);
        System.out.println("\nPrevious pathway FIs: " + previousPathwayFIs.size());
        Set<String> previousPredictedFIs = fu.loadInteractions(dirName + "FIs_Predicted_043009.txt");
        previousPredictedFIs = replaceSpaceWithTab(previousPredictedFIs);
        System.out.println("Previous predicted FIs: " + previousPredictedFIs.size());
        Set<String> previousFIs = new HashSet<String>(previousPathwayFIs);
        previousFIs.addAll(previousPredictedFIs);
        System.out.println("Total: " + previousFIs.size());
        // Check if previous pathway FIs have been contained by current FIs
        Set<String> sharedFIs = InteractionUtilities.getShared(currentPathwayFIs, previousPathwayFIs);
        double percent = (double) sharedFIs.size() / Math.min(currentPathwayFIs.size(), previousPathwayFIs.size());
        System.out.println("\nShared pathway FIs: " + sharedFIs.size() + " (" + percent + ")");
        sharedFIs = InteractionUtilities.getShared(currentPredictedFIs, previousPredictedFIs);
        percent = (double) sharedFIs.size() / Math.min(currentPredictedFIs.size(), previousPredictedFIs.size());
        System.out.println("Shared predicted FIs: " + sharedFIs.size() + " (" + percent + ")");
        sharedFIs = InteractionUtilities.getShared(currentFIs, previousFIs);
        percent = (double) sharedFIs.size() / Math.min(currentFIs.size(), previousFIs.size());
        System.out.println("Shared total FIs: " + sharedFIs.size() + " (" + percent + ")");
        
        // Check sharing for predictedFIs
        checkFeaturesInPathwayFIs(currentPathwayFIs,
                                  previousPathwayFIs,
                                  previousPredictedFIs,
                                  "Previous predicted FIs");
        
        Set<String> carlosGeneExp = fu.loadInteractions(dirName + "CarlosCoExp_Norm.txt");
        carlosGeneExp = replaceSpaceWithTab(carlosGeneExp);
        checkFeaturesInPathwayFIs(currentPathwayFIs, 
                                  previousPathwayFIs,
                                  carlosGeneExp, 
                                  "Carlos Gene Expression");
        
        Set<String> pavlidisGeneExp = fu.loadInteractions(dirName + "PavlidisCoExp_Norm.txt");
        pavlidisGeneExp = replaceSpaceWithTab(pavlidisGeneExp);
        checkFeaturesInPathwayFIs(currentPathwayFIs, 
                                  previousPathwayFIs, 
                                  pavlidisGeneExp, 
                                  "Pavlidis Gene Expression");
        
        // Check for old version of human PPIs
        Set<String> humanPPIs = new HashSet<String>();
        String fileName = dirName + "HumanPPIs_Less4_intact.txt";
        humanPPIs.addAll(fu.loadInteractions(fileName));
        fileName = dirName + "HumanPPIs_BioGrid.txt";
        humanPPIs.addAll(fu.loadInteractions(fileName));
        fileName = dirName + "HumanPPIs_HPRD.txt";
        humanPPIs.addAll(fu.loadInteractions(fileName));
        humanPPIs = replaceSpaceWithTab(humanPPIs);
        checkFeaturesInPathwayFIs(currentPathwayFIs, 
                                  previousPathwayFIs,
                                  humanPPIs,
                                  "Previous Human PPIs");
        
        Set<String> sharedHumanPPIsAndPrePredictedFIs = InteractionUtilities.getShared(humanPPIs, previousFIs);
        checkFeaturesInPathwayFIs(currentPathwayFIs, 
                                  previousPathwayFIs, 
                                  sharedHumanPPIsAndPrePredictedFIs, 
                                  "Shared Human PPIs and Previous Predicted FIs");
        
        Set<String> flyPPIs = fu.loadInteractions(dirName + "humanPPIsFromFlyInUniProt_Norm.txt");
        flyPPIs = replaceSpaceWithTab(flyPPIs);
        checkFeaturesInPathwayFIs(currentPathwayFIs,
                                  previousPathwayFIs,
                                  flyPPIs,
                                  "Previous Fly PPIs");

        Set<String> wormPPIs = fu.loadInteractions(dirName + "humanPPIsFromWormInUniProt_Norm.txt");
        wormPPIs = replaceSpaceWithTab(wormPPIs);
        checkFeaturesInPathwayFIs(currentPathwayFIs,
                                  previousPathwayFIs,
                                  wormPPIs,
                                  "Previous Worm PPIs");
        
        Set<String> yeastPPIs = fu.loadInteractions(dirName + "humanPPIsFromYeastInUniProt_Norm.txt");
        yeastPPIs = replaceSpaceWithTab(yeastPPIs);
        checkFeaturesInPathwayFIs(currentPathwayFIs,
                                  previousPathwayFIs,
                                  yeastPPIs,
                                  "Previous Yeast PPIs");
        
        // Check for PPIs using the current version: actually it will be better to use the old version
        humanPPIs = fu.loadInteractions(FIConfiguration.getConfiguration().get("IREFINDEX_HUMAN_PPI_FILE"));
        checkFeaturesInPathwayFIs(currentPathwayFIs, 
                                  previousPathwayFIs,
                                  humanPPIs, 
                                  "Current human PPIs");
        
//        previousPredictedFIs.removeAll(sharedFIs);
//        System.out.println("\nNot shared predicted FIs in previous:");
//        Iterator<String> it = previousPredictedFIs.iterator(); 
//        
//        // Check with current NBC
//        NBCAnalyzer nbcAnalyzer = new NBCAnalyzer();
//        NaiveBayesClassifier classifier = nbcAnalyzer.loadSavedNBC();
//        // Get pairs from all features
//        Set<String> allPairs = nbcAnalyzer.loadPairForPrediction();
//        FeatureHandlerForV3 featureHandler = new FeatureHandlerForV3();
//        Map<String, PositiveChecker> featureToChecker = featureHandler.loadFeatureToChecker();
//        for (int i = 0; i < 100; i++) {
//            String fi = it.next();
//            System.out.println(fi);
//            Double score = classifier.calculateScore(fi, featureToChecker);
//            System.out.println("Score: " + score);
//        }
    }

    public void checkFeaturesInPathwayFIs(final Set<String> currentPathwayFIs,
                                          Set<String> previousPathwayFIs,
                                          Set<String> humanPPIs,
                                          String featureName) {
        Set<String> sharedFIs;
        double percent;
        System.out.println("\nChecking " + featureName);
        System.out.println("Total pairs: " + humanPPIs.size());
        humanPPIs.removeAll(previousPathwayFIs);
        System.out.println("After removing old pathwayFIs: " + humanPPIs.size());
        sharedFIs = InteractionUtilities.getShared(humanPPIs, currentPathwayFIs);
        percent = (double) sharedFIs.size() / humanPPIs.size();
        System.out.println("Pairs in current pathway FIs: " + sharedFIs.size() + " (" + percent + ")");
        
        FeatureChecker featureChecker = new FeatureChecker();
        PositiveChecker positiveChecker = new PositiveChecker() {
            
            @Override
            public boolean isPositive(String pair) {
                return currentPathwayFIs.contains(pair);
            }
        };
        featureChecker.checkFeatureOddsRatio(humanPPIs, positiveChecker);
    }
    
    private Set<String> replaceSpaceWithTab(Set<String> fis) {
        Set<String> rtn = new HashSet<String>();
        for (String fi : fis) {
            int index = fi.indexOf(" ");
            if (index > 0)
                rtn.add(fi.substring(0, index) + "\t" + fi.substring(index + 1));
            else
                rtn.add(fi);
        }
        return rtn;
    }
}
