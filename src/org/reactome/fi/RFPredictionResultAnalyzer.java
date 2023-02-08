package org.reactome.fi;

import java.io.File;
import java.io.IOException;
import java.lang.reflect.Method;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.apache.log4j.Logger;
import org.gk.persistence.MySQLAdaptor;
import org.junit.Test;
import org.reactome.data.ReactomeAnalyzer;
import org.reactome.data.UniProtAnalyzer;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.InteractionUtilities;
import org.reactome.funcInt.Evidence;
import org.reactome.r3.util.FileUtility;

import com.google.common.base.Functions;

import tech.tablesaw.api.FloatColumn;
import tech.tablesaw.api.IntColumn;
import tech.tablesaw.api.Row;
import tech.tablesaw.api.Table;
import tech.tablesaw.plotly.Plot;
import tech.tablesaw.plotly.api.LinePlot;

/**
 * This class is used to process and analyze the prediction results from the trained RF coming from
 * the fi_network_ml project: https://github.com/reactome-idg/fi-network-ml. 
 * @author wug
 *
 */
public class RFPredictionResultAnalyzer {
    private final static Logger logger = Logger.getLogger(RFPredictionResultAnalyzer.class);
    
    public RFPredictionResultAnalyzer() {
    }
    
    /**
     * Use this method to generate the three FI files in genes:
     * FIsInGene, FIsInGenes_Pathways, FIsInGene_Predicted.
     * @throws Exception
     */
    @Test
    public void generateFIFiles() throws Exception {
        double threshold = Double.parseDouble(FIConfiguration.getConfiguration().get("CUT_OFF_VALUE"));
        Map<String, Double> fi2score = loadPredictedFIsWithScores();
        Set<String> predictedFIs = new HashSet<>();
        for (String fi : fi2score.keySet()) {
            Double score = fi2score.get(fi);
            if (score >= threshold)
                predictedFIs.add(fi);
        }
        reportNumbers(predictedFIs, "Predicted FIs in Genes for " + threshold, false);
        predictedFIs = normalizeGeneNames(predictedFIs);
        reportNumbers(predictedFIs, "After name normalization", false);
        // Get the annotated FIs
        Set<String> annotatedFIs = new FIFileAnalyzer().loadPathwayAndTFTargetFIs();
        reportNumbers(annotatedFIs, "Annotated FIs in UniProt", true);
        // Convert these FIs from UniProt to genes
        Set<String> annotatedFIsInGenes = convertFIsToGenes(annotatedFIs);
        reportNumbers(annotatedFIsInGenes, "Annotated FIs in Genes", false);
        // There might be some inconsistence for FIs collected from other resources.
        Set<String> annotatedFIsInGenes1 = normalizeGeneNames(annotatedFIsInGenes);
        reportNumbers(annotatedFIsInGenes1, "After name normalization", false);
        // Just a check
        if (annotatedFIsInGenes.size() != annotatedFIsInGenes1.size())
            throw new IllegalStateException("Annotated FIs have used not normalized names!");
        
        // Need to remove annotated FIs from normalized predicted FIs
        // Annotated FIs should be normalized already
        predictedFIs.removeAll(annotatedFIsInGenes); 
        reportNumbers(predictedFIs, "Remove annotateFIs in normalized predicted FIs", false);
        
        Set<String> totalFIs = new HashSet<>(predictedFIs);
        totalFIs.addAll(annotatedFIsInGenes);
        reportNumbers(totalFIs, "Total FIs in Genes", false);
        // Output these files
        logger.info("Saving FIs...");
        String fileName = FIConfiguration.getConfiguration().get("GENE_FI_FILE_NAME");
        FileUtility fu = new FileUtility();
        fu.saveInteractions(totalFIs, fileName);
        fileName = FIConfiguration.getConfiguration().get("GENE_FI_PATHWAY_FILE_NAME");
        fu.saveInteractions(annotatedFIs, fileName);
        fileName = FIConfiguration.getConfiguration().get("GENE_FI_PREDICTED_FILE_NAME");
        fu.saveInteractions(predictedFIs, fileName);
        logger.info("Done.");
    }
    
    /**
     * Make sure all gene names are current based on the current Reactome release.
     * @param fis
     * @return
     * @throws Exception
     */
    private Set<String> normalizeGeneNames(Set<String> fis) throws Exception {
        Map<String, String> synonym2gene = getSynonym2Gene();
        Set<String> rtn = new HashSet<>();
        for (String fi : fis) {
            String fi1 = normalizeFI(fi, synonym2gene);
            if (fi1 == null)
                continue;
            rtn.add(fi1);
        }
        return rtn;
    }

    private String normalizeFI(String fi, Map<String, String> synonym2gene) {
        String[] genes = fi.split("\t");
        String gene1 = synonym2gene.get(genes[0]);
        String gene2 = synonym2gene.get(genes[1]);
        // Don't want to include self interaction
        if (gene1.equals(gene2))
            return null;
        String fi1 = InteractionUtilities.generateFIFromGene(gene1, gene2);
        return fi1;
    }
    
    @Test
    public void countLastYearFINumbers() throws Exception {
        String resultDir = FIConfiguration.getConfiguration().get("RESULT_DIR");
        resultDir = resultDir + "/../../fi_network_build/2021/";
        File dir = new File(resultDir);
        if (!dir.exists())
            throw new IllegalStateException("Cannot find: " + dir.getAbsolutePath());
        Set<String> pathwayFIs = new HashSet<>();
        FileUtility fu = new FileUtility();
        for (File file : dir.listFiles()) {
            String fileName = file.getName();
            if (fileName.startsWith("FIs_")) {
                logger.info("Loading " + fileName + "...");
                Set<String> fis = fu.loadInteractions(file.getAbsolutePath());
                logger.info("FIs: " + fis.size());
                pathwayFIs.addAll(fis);
            }
        }
        reportNumbers(pathwayFIs, "Pathway FIs", true);
        String fileName = resultDir + "PredictedFIs_122921.txt";
        Set<String> predictedFIs = fu.loadInteractions(fileName);
        reportNumbers(predictedFIs, "Predicted FIs", true);
        Set<String> totalFIs = new HashSet<>(pathwayFIs);
        totalFIs.addAll(predictedFIs);
        reportNumbers(totalFIs, "Total FIs", true);
        // Numbers in genes
        fileName = resultDir + "FIsInGene_Pathway_122921.txt";
        Set<String> pathwayFIsInGenes = fu.loadInteractions(fileName);
        reportNumbers(pathwayFIsInGenes, "Pathway FIs in Genes", false);
        fileName = resultDir + "FIsInGene_Predicted_122921.txt";
        Set<String> predictedFIsInGenes = fu.loadInteractions(fileName);
        reportNumbers(predictedFIsInGenes, "Predicted FIs in Genes", false);
        fileName = resultDir + "FIsInGene_122921.txt";
        Set<String> totalFIsInGenes = fu.loadInteractions(fileName);
        reportNumbers(totalFIsInGenes, "Total FIs in Genes", false);
    }
    
    @Test
    public void checkAnnotatedFIs() throws Exception {
        // Get the annotated FIs
        Set<String> annotatedFIs = new FIFileAnalyzer().loadPathwayAndTFTargetFIs();
        reportNumbers(annotatedFIs, "Annotated FIs", true);
        // Convert these FIs from UniProt to genes
        Set<String> annotatedFIsInGenes = convertFIsToGenes(annotatedFIs);
        reportNumbers(annotatedFIsInGenes, "Annotated FIs in Genes", false);
    }
    
    private void reportNumbers(Set<String> fis, 
                               String title,
                               boolean filterBySwissProtIds) throws IOException {
        StringBuilder builder = new StringBuilder();
        builder.append(title).append(": ");
        builder.append("\n\tFI: " + fis.size());
        Set<String> ids = InteractionUtilities.grepIDsFromInteractions(fis);
        if (filterBySwissProtIds) {
            Set<String> swissProtIds = new UniProtAnalyzer().loadSwissProtIds();
            logger.info("Total SwissProt ids: " + swissProtIds.size());
            ids.retainAll(swissProtIds);
        }
        builder.append("\n\tIds: " + ids.size());
        logger.info(builder.toString());
    }
    
    @Test
    public void checkThreshold() throws Exception {
//        Starting with this version, the FI network will be based on genes since they are more likely used
//        in many data analysis and visualization even though the network is actually protein-based.
        // The prediction file is in gene names. However, to follow the old FI build procedures, we need to 
        // use UniProt ids. 
        Map<String, Double> fi2score = loadPredictedFIsWithScores();
        plotPredictedNumberVsThreshold(fi2score);
    }

    private Map<String, Double> loadPredictedFIsWithScores() throws Exception {
        String fileName = FIConfiguration.getConfiguration().get("RF_PREDICTION_FILE");
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        Map<String, Double> fi2score = new HashMap<>();
        String line = fu.readLine(); // The header
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split(",");
            // Seven columns: GenePair,rep1,rep2,rep3,mean,stdev,%cv. %cv = 100 * mean/stdev
            // Need to reformat the FIs
            String[] geneTokens = tokens[0].split("_");
            String fiInGenes = InteractionUtilities.generateFIFromGene(geneTokens[0], 
                                                                       geneTokens[1]);
            fi2score.put(fiInGenes, Double.parseDouble(tokens[4])); // Use mean for our analysis
        }
        fu.close();
        logger.info("Total predicted pairs: " + fi2score.size());
        // Get the annotated FIs
        Set<String> annotatedFIs = new FIFileAnalyzer().loadPathwayAndTFTargetFIs();
        // Convert these FIs from UniProt to genes
        Set<String> annotatedFIsInGenes = convertFIsToGenes(annotatedFIs);
        // Remove these UniProt FIs for quick calculation below
        fi2score.keySet().removeAll(annotatedFIsInGenes);
        logger.info("Total predicted pairs after removing annotated FIs: " + fi2score.size());
        return fi2score;
    }
    
    private void plotPredictedNumberVsThreshold(Map<String, Double> fi2score) throws IOException {
        Table table = Table.create();
        String thresholdColName = "Score";
        FloatColumn thresholdCol = FloatColumn.create(thresholdColName);
        String numberColName = "Number of Predicted FIs";
        IntColumn numberCol = IntColumn.create(numberColName);
        table.addColumns(thresholdCol, numberCol);
        String numberOfGeneColName = "Number of Genes";
        IntColumn geneNumberCol = IntColumn.create(numberOfGeneColName);
        table.addColumns(geneNumberCol);
        // Create 100 points
        for (float threshold = 0.0f; threshold <= 1.0f; threshold += 0.01) {
            logger.info("Threshold " + threshold + "...");
            Row row = table.appendRow();
            row.setFloat(thresholdColName, threshold);
            int[] numbers = countPredictedFIsAndGenes(fi2score, threshold);
            row.setInt(numberColName, numbers[0]);
            row.setInt(numberOfGeneColName, numbers[1]);
        }
        logger.info("Finished table building.");
        // Output for better plot
        String plotFileName = FIConfiguration.getConfiguration().get("RESULT_DIR") + 
                              File.separator + 
                              "NumberOfPredictedFIsVsThreshold_02042022.csv";
        table.write().csv(plotFileName);
        logger.info("Table saved.");
        Plot.show(LinePlot.create("Number of Predicted FIs vs Score",
                table,
                thresholdColName, 
                numberColName));
        Plot.show(LinePlot.create("Number of Genes vs Score",
                table,
                thresholdColName,
                numberOfGeneColName));
        logger.info("All done.");
    }
    
    /**
     * Count predicted FIs and genes.
     * @param fi2score
     * @param threshold
     * @return two elements, the first is number of  Fis and second number of genes
     */
    private int[] countPredictedFIsAndGenes(Map<String, Double> fi2score,
                                            float threshold) {
        int counter = 0;
        Set<String> genes = new HashSet<>();
        int index = 0;
        for (String fi : fi2score.keySet()) {
            Double score = fi2score.get(fi);
            if (score >= threshold) {
                counter ++;
                // Use index should be faster than split.
                index = fi.indexOf('\t');
                genes.add(fi.substring(0, index));
                genes.add(fi.substring(index + 1));
            }
        }
        return new int[]{counter, genes.size()};
    }
    
    private Set<String> convertFIsToGenes(Set<String> fis) throws Exception {
        Map<String, String> uniprot2gene = getUniprot2Gene();
        Set<String> fisInGenes = new HashSet<>();
        for (String fi : fis) {
            String[] tokens = fi.split("\t");
            String gene1 = uniprot2gene.get(tokens[0]);
            String gene2 = uniprot2gene.get(tokens[1]);
            // Don't want to keep the self interaction
            if (gene1 == null || gene2 == null || gene1.equals(gene2))
                continue;
            String fiInGene = InteractionUtilities.generateFIFromGene(gene1, gene2);
            fisInGenes.add(fiInGene);
        }
        return fisInGenes;
    }
    
    private Map<String, String> getUniprot2Gene() throws Exception {
        ReactomeAnalyzer analyzer = new ReactomeAnalyzer();
        MySQLAdaptor dba = (MySQLAdaptor) analyzer.getMySQLAdaptor();
        Map<String, String> uni2gene = analyzer.getUniProtToGeneMap(dba);
        logger.info("Size of uni2gene: " + uni2gene.size());
        return uni2gene;
    }
    
    public Map<String, String> getGene2UniProt() throws Exception {
        Map<String, String> uniprot2gene = getUniprot2Gene();
        Map<String, Set<String>> gene2uniprot = new HashMap<>();
        uniprot2gene.forEach((u, g) -> {
            if (g == null)
                return; // Do nothing if no gene there
            Set<String> set = gene2uniprot.get(g);
            if (set == null) {
                set = new HashSet<>();
                gene2uniprot.put(g, set);
            }
            set.add(u);
        });
        Map<String, String> rtnGene2Uniprot = new HashMap<>();
        gene2uniprot.forEach((g, set) -> {
            // Choose the first accession number that are not there alphabetically
            List<String> list = set.stream().sorted().collect(Collectors.toList());
            rtnGene2Uniprot.put(g, list.get(0));
        });
        // Make sure the mapping is unique
        Set<String> valueSet = new HashSet<>(rtnGene2Uniprot.values());
        if (rtnGene2Uniprot.keySet().size() != valueSet.size())
            throw new IllegalStateException("Cannot get one to one mapping between gene and uniprot!");
        logger.info("The size of gene2uniprot map: " + rtnGene2Uniprot.size());
        return rtnGene2Uniprot;
    }
    
    @Test
    public void testGetGene2UniProt() throws Exception {
        Map<String, String> gene2uni = getGene2UniProt();
        // Check if all predicted Fis in genes can be mapped
        String fileName = FIConfiguration.getConfiguration().get("GENE_FI_PREDICTED_FILE_NAME");
        Set<String> predictedFIs = new FileUtility().loadInteractions(fileName);
        System.out.println("Predicted FI: " + predictedFIs.size());
        Set<String> FIs = new HashSet<>();
        for (String fiInGenes : predictedFIs) {
            String[] genes = fiInGenes.split("\t");
            String acc1 = gene2uni.get(genes[0]);
            String acc2 = gene2uni.get(genes[1]);
            if (acc1 == null || acc2 == null) {
                System.out.println(fiInGenes + " cannot be mapped to UniProt!");
                continue;
            }
            String fi = InteractionUtilities.generateFIFromGene(acc1, acc2);
            FIs.add(fi);
        }
        System.out.println("Total FIs in UniProt: " + FIs.size());
    }
    
    public Map<String, String> mapFIsInGenesToFIsInUniProt(Set<String> fisInGenes) throws Exception {
        Map<String, String> gene2uni = getGene2UniProt();
        Map<String, String> map = new HashMap<>();
        for (String fiInGenes : fisInGenes) {
            String[] genes = fiInGenes.split("\t");
            String acc1 = gene2uni.get(genes[0]);
            String acc2 = gene2uni.get(genes[1]);
            if (acc1 == null || acc2 == null) {
                logger.error(fiInGenes + " cannot be mapped to UniProt!");
                continue;
            }
            String fi = InteractionUtilities.generateFIFromGene(acc1, acc2);
            map.put(fiInGenes, fi);
        }
        Set<String> fisInUniProt = new HashSet<>(map.values());
        if (fisInGenes.size() != fisInUniProt.size())
            throw new IllegalStateException("Cannot map FIs in genes to FIs in UniProt one to one!");
        return map;
    }
    
    public Map<String, String> getSynonym2Gene() throws Exception {
        ReactomeAnalyzer analyzer = new ReactomeAnalyzer();
        MySQLAdaptor dba = (MySQLAdaptor) analyzer.getMySQLAdaptor();
        Map<String, String> synonym2gene = analyzer.getSynonymToGeneName(dba);
        logger.info("Size of synonym2gene: " + synonym2gene.size());
        return synonym2gene;
    }
    
    @Test
    public void testGetUniprot2Gene() throws Exception {
        Map<String, String> uni2gene = getUniprot2Gene();
        Set<String> genes = uni2gene.values().stream().filter(g -> g != null).collect(Collectors.toSet());
        logger.info("Total genes: " + genes.size());
        Set<String> noGenesUniProt = uni2gene.keySet().stream().filter(u -> uni2gene.get(u) == null).collect(Collectors.toSet());
        double ratio = (double)noGenesUniProt.size() / uni2gene.size();
        logger.info("Total UniProt ids having no genes: " + noGenesUniProt.size() + " (" + ratio +  ")");
        noGenesUniProt.forEach(System.out::println);
    }
    
    @Test
    public void checkGeneToUniProtMapInPredictedFIs() throws Exception {
        String fileName = FIConfiguration.getConfiguration().get("GENE_FI_PREDICTED_FILE_NAME");
        Set<String> predictedFIs = new FileUtility().loadInteractions(fileName);
        logger.info("Total predicted fis: " + predictedFIs.size());
        Map<String, String> uniprot2gene = getUniprot2Gene();
        Map<String, Set<String>> gene2uniprot = new HashMap<>();
        for (String uniprot : uniprot2gene.keySet()) {
            String gene = uniprot2gene.get(uniprot);
            gene2uniprot.compute(gene, (key, set) -> {
               if (set == null)
                   set = new HashSet<>();
               set.add(uniprot);
               return set;
            });
        }
        Set<String> genesWithMoreUniProts = new HashSet<>();
        Set<String> genesWithoutUniProts = new HashSet<>();
        for (String fi : predictedFIs) {
            String[] genes = fi.split("\t");
            Set<String> uniprots = gene2uniprot.get(genes[0]);
            if (uniprots == null)
                genesWithoutUniProts.add(genes[0]);
            else if (uniprots.size() > 1)
                genesWithMoreUniProts.add(genes[0]);
            uniprots = gene2uniprot.get(genes[1]);
            if (uniprots == null)
                genesWithoutUniProts.add(genes[1]);
            else if (uniprots.size() > 1)
                genesWithMoreUniProts.add(genes[1]);
        }
        logger.info("Genes having more than one UniProts: " + genesWithMoreUniProts.size());
        for (String gene : genesWithMoreUniProts) {
            logger.info(gene + ": " + gene2uniprot.get(gene));
        }
        logger.info("Genes having no UniProts: " + genesWithoutUniProts.size());
        genesWithoutUniProts.forEach(System.out::println);
    }
    
    /**
     * The following method is used to check the feature file to generate some related Java 
     * properties and meta information for ReactomeFIViz.
     * @throws IOException
     */
    @Test
    public void checkFeatureFileHeaders() throws IOException {
        String featureFile = FIConfiguration.getConfiguration().get("RF_FEATURE_FILE");
        FileUtility fu = new FileUtility();
        fu.setInput(featureFile);
        String header = fu.readLine();
        fu.close();
        String[] tokens = header.split(",");
        List<String> features = Stream.of(tokens).skip(1).collect(Collectors.toList());
        logger.info("Total features: " + features.size() + "\n");
        Function<String, String> format = feature -> {
            return getJavaPropertyName(feature);
        };
        // For Java Evidence class
        features.stream().forEach(feature -> {
            feature = format.apply(feature);
            logger.info("private Boolean " + feature + ";");
        });
        // For Hibernate mapping, which uses the old version with XML. 
        logger.info("\nHiberante Mapping XML:\n");
        features.stream().forEach(feature -> {
            feature = format.apply(feature);
            logger.info("<property name=\"" + feature + "\" />");
        });
    }

    private String getJavaPropertyName(String feature) {
        if (!feature.startsWith("GTE") && !feature.startsWith("TCGA"))
            feature = feature.substring(0, 1).toLowerCase() + feature.substring(1);
        feature = feature.replaceAll("-", "_");
        return feature;
    }
    
    /**
     * Load the predicted FIs together with evidence. This method should be used to dump predicted FIs into
     * a database together with Evidence. The other method, generateFIFiles(), should be used to generate
     * the FIs files without evidence. The returned FIs from this method have been normalized for gene names.
     * @param fis
     * @return
     * @throws IOException
     */
    public Map<String, Evidence> loadPredictedFIsWithEvidence() throws Exception {
        // Load the scores
        Map<String, Double> fi2score = loadPredictedFIsWithScores();
        // Do filtering
        double threshold = Double.parseDouble(FIConfiguration.getConfiguration().get("CUT_OFF_VALUE"));
        Map<String, Double> filteredFI2score = fi2score.keySet().stream()
                .filter(fi -> (fi2score.get(fi) >= threshold))
                .collect(Collectors.toMap(Functions.identity(), fi -> fi2score.get(fi)));
        logger.info("Total filtered fis2score: " + filteredFI2score.size());
        
        String featureFile = FIConfiguration.getConfiguration().get("RF_FEATURE_FILE");
        FileUtility fu = new FileUtility();
        fu.setInput(featureFile);
        String header = fu.readLine();
        String[] features = header.split(",");
        String line = null;
        String[] tokens = null;
        int index = 0;
        String fiInFeature = null;
        Map<String, Evidence> fi2evidence = new HashMap<>();
        while ((line = fu.readLine()) != null) {
            // Quick look at the FI in the file first
            index = line.indexOf(",");
            fiInFeature = line.substring(0, index);
            if (!filteredFI2score.keySet().contains(fiInFeature))
                continue;
            // A slower step
            tokens = line.split(",");
            Evidence evidence = new Evidence();
            evidence.setScore(filteredFI2score.get(tokens[0]));
            for (int i = 1; i < tokens.length; i++) {
                // Use Java reflection to assign the feature
                // Check method 
                String propName = getJavaPropertyName(features[i]);
                String methodName = null;
                if (Character.isUpperCase(propName.charAt(1)))
                    methodName = "set" + propName; // Use the propName directly (e.g. gOBPSharing, TCGAXX)
                else
                    methodName = "set" + propName.substring(0, 1).toUpperCase() + propName.substring(1);
                Method method = Evidence.class.getMethod(methodName, Boolean.class);
                method.invoke(evidence, tokens[i].equals("1") ? Boolean.TRUE : Boolean.FALSE);
            }
            fi2evidence.put(tokens[0], evidence);
        }
        fu.close();
        logger.info("The size of fi2evidence: " + fi2evidence.size());
        // Get the annotated FIs so that we can remove normalized predicted FIs that may be matched to annotated FIs
        Set<String> annotatedFIs = new FIFileAnalyzer().loadPathwayAndTFTargetFIs();
        // Convert these FIs from UniProt to genes
        Set<String> annotatedFIsInGenes = convertFIsToGenes(annotatedFIs);
        logger.info("Annotated FIs in genes: " + annotatedFIsInGenes.size());
        // Need to normalize
        Map<String, String> synonym2gene = getSynonym2Gene();
        Map<String, Evidence> rtnFi2evidence = new HashMap<>();
        for (String fi : fi2evidence.keySet()) {
            Evidence evidence = fi2evidence.get(fi);
            String normalizedFI = normalizeFI(fi, synonym2gene);
            if (normalizedFI == null || annotatedFIsInGenes.contains(normalizedFI))
                continue; // Remove self interaction or annotated FIs
            // Some of genes may be mapped to the other FIs because
            // of this normalization. Just use whatever evidence there.
            if (rtnFi2evidence.containsKey(normalizedFI)) {
                logger.info("Warning: " + normalizedFI + " is in the rtnFi2evidence map already!");
                continue;
            }
            rtnFi2evidence.put(normalizedFI, evidence);
        }
        logger.info("After name normalization: " + rtnFi2evidence.size());
        return rtnFi2evidence;
    }
    
    @Test
    public void testLoadPredictedFIsWithEvidence() throws Exception {
        Map<String, Evidence> fi2evidence = loadPredictedFIsWithEvidence();
        String fileName = FIConfiguration.getConfiguration().get("GENE_FI_PREDICTED_FILE_NAME");
        Set<String> predictedFIs = new FileUtility().loadInteractions(fileName);
        System.out.println("Predicted FIs in file: " + predictedFIs.size());
        Set<String> fisInEvidence = new HashSet<>(fi2evidence.keySet());
        System.out.println("Predicted FIs with evidence: " + fisInEvidence.size());
        fisInEvidence.removeAll(predictedFIs);
        System.out.println("Removed all FIs in file: " + fisInEvidence.size());
        fisInEvidence.forEach(System.out::println);
    }
    
    @Test
    public void testNormalizeGeneNames() throws Exception {
        Set<String> fis = new HashSet<>();
        fis.add("IARS1\tQARS1");
        Set<String> normalized = normalizeGeneNames(fis);
        System.out.println(normalized);
    }
    
}
