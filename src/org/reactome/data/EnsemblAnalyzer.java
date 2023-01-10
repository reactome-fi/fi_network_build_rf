/*
 * Created on Nov 28, 2006
 *
 */
package org.reactome.data;

import java.io.IOException;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;
import org.junit.Test;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.FileUtility;

/**
 * This class is used to query ensembl MySQL database.
 * @author guanming
 *
 */
public class EnsemblAnalyzer {
    private static final Logger logger = Logger.getLogger(EnsemblAnalyzer.class);
    
    private Connection getConnection() throws Exception {
        String dbName = FIConfiguration.getConfiguration().get("ENSEMBL_COMPARA_DATABASE");
        return FIConfiguration.getConnection(dbName);
    }
    
    /**
     * This method is used to dump protein families for four organisms.
     * @throws Exception
     */
    @Test
    public void dumpProteinFamilies() throws Exception {
        // Pick uniprot identifiers only
        // As of December, 2014, table member has been replaced by seq_member
        String queryString = "select f.family_id, m.stable_id from family_member f," +
        		" seq_member m where f.seq_member_id = m.seq_member_id and m.taxon_id = ? " +
        		"and m.source_name like 'UniProt%'";
        Connection connection = getConnection();
        connection.setReadOnly(true);
        PreparedStatement stat = connection.prepareStatement(queryString);
        List<Integer> neededTaxonIds = getNeededTaxonIds();
        Map<String, Set<String>> familyToProteins = new HashMap<String, Set<String>>();
        for (Integer taxonId : neededTaxonIds) {
            logger.info("Dump " + taxonId + "...");
            stat.setInt(1, taxonId);
            ResultSet resultset = stat.executeQuery();
//            int count = 0;
            while (resultset.next()) {
                String family = resultset.getString(1);
                String uniProtId = resultset.getString(2);
                Set<String> set = familyToProteins.get(family);
                if (set == null) {
                    set = new HashSet<String>();
                    familyToProteins.put(family, set);
                }
                set.add(taxonId + ":" + uniProtId);
//                System.out.println("count: " + count++);
            }
            resultset.close();
            logger.info("Finish taxon: " + taxonId);
//            System.out.println("Finish taxon: " + taxonId);
        }
        stat.close();
        connection.close();
        logger.info("Size of familyToProteins(): " + familyToProteins.size());
        filterFamilies(familyToProteins);
        FileUtility fu = new FileUtility();
        fu.saveSetMapInSort(familyToProteins, 
                            FIConfiguration.getConfiguration().get("ENSEMBL_PROTEIN_FAMILIES"));
    }
    
    public Map<String, Set<String>> loadYeastToHumanMapInUniProt() throws IOException {
        return loadToHumanMapInUniProt("559292"); // need to use this yeast type!!!
//        return loadToHumanMapInUniProt("4932");
    }
    
    public Map<String, Set<String>> loadWormToHumanMapInUniProt() throws IOException {
        return loadToHumanMapInUniProt("6239");
    }
    
    public Map<String, Set<String>> loadFlyToHumanMapInUniProt() throws IOException {
        return loadToHumanMapInUniProt("7227");
    }
    
    public Map<String, Set<String>> loadMouseToHumanMapInUniProt() throws IOException {
        return loadToHumanMapInUniProt("10090");
    }
    
    private Map<String, Set<String>> loadToHumanMapInUniProt(String taxonId) throws IOException {
        FileUtility fu = new FileUtility();
        String fileName = FIConfiguration.getConfiguration().get("ENSEMBL_PROTEIN_FAMILIES");
        Map<String, Set<String>> familyToProteins = fu.loadSetMap(fileName);
        Set<String> humanIds = new HashSet<String>();
        Set<String> otherIds = new HashSet<String>();
        // To be returned
        Map<String, Set<String>> map = new HashMap<String, Set<String>>();
        for (String family : familyToProteins.keySet()) {
            Set<String> proteins = familyToProteins.get(family);
            humanIds.clear();
            otherIds.clear();
            splitIds(proteins, humanIds, otherIds, taxonId);
            for (String otherId : otherIds) {
                Set<String> humanSet = map.get(otherId);
                if (humanSet == null) {
                    humanSet = new HashSet<String>();
                    map.put(otherId, humanSet);
                }
                humanSet.addAll(humanIds);
            }
        }
        return map;
    }
    
    private void splitIds(Set<String> proteins,
                          Set<String> humanIds,
                          Set<String> otherIds,
                          String taxonId) {
        for (String protein : proteins) {
            if (protein.startsWith("9606:")) {
                // This is a human protein
                humanIds.add(protein.substring(5));
            }
            else if (protein.startsWith(taxonId)) {
                otherIds.add(protein.substring(taxonId.length() + 1)); // 1 for ":".
            }
        }
    }
     
    
    /**
     * Filter protein families so that a family containing at least two species and one of them
     * should be homo sapiens.
     * @param familyToProteins
     */
    private void filterFamilies(Map<String, Set<String>> familyToProteins) {
        System.out.println("Total families before filtering: " + familyToProteins.size());
        for (Iterator<String> it = familyToProteins.keySet().iterator(); it.hasNext();) {
            String family = it.next();
            Set<String> proteins = familyToProteins.get(family);
            Set<String> species = extractSpecies(proteins);
            if (species.size() == 1 ||
                !species.contains("9606")) {
                it.remove();
            }
        }
        System.out.println("Total families after filtering: " + familyToProteins.size());
    }
    
    private Set<String> extractSpecies(Set<String> proteins) {
        Set<String> species = new HashSet<String>();
        int index = 0;
        for (String protein : proteins) {
            index = protein.indexOf(":");
            species.add(protein.substring(0, index));
        }
        return species;
    }
    
    /**
     * This method is used to filter a big sequence.txt file downloaded from ensembl
     * compara to sequence files used for four species only.
     * @throws Exception
     */
    public void filterSequences() throws Exception {
        List<Long> sequenceIds = getSequenceIdsForNeededTaxons();
        Set<Long> idSet = new HashSet<Long>(sequenceIds); // Should be fast
        System.out.println("Total ids: " + idSet.size());
        String inFile = FIConfiguration.getConfiguration().get("ENSEMBL_DIR") + "sequence.txt";
        String outFile = FIConfiguration.getConfiguration().get("ENSEMBL_DIR") + "sequence_filtered.txt";
        FileUtility inFu = new FileUtility();
        inFu.setInput(inFile);
        FileUtility outFu = new FileUtility();
        outFu.setOutput(outFile);
        String line = null;
        int index = 0;
        String id = null;
        while ((line = inFu.readLine()) != null) {
            index = line.indexOf("\t");
            id = line.substring(0, index);
            if(idSet.contains(Long.parseLong(id)))
                outFu.printLine(line);
        }
        inFu.close();
        outFu.close();
    }
    
    private List<Long> getSequenceIdsForNeededTaxons() throws Exception {
        Connection connection = getConnection();
        List<Integer> taxonIds = getNeededTaxonIds();
        List<Long> sequenceIds = new ArrayList<Long>();
        String query = "SELECT sequence_id FROM member WHERE taxon_id = ?";
        PreparedStatement stat = connection.prepareStatement(query);
        for (Integer id : taxonIds) {
            stat.setInt(1, id);
            ResultSet resultSet = stat.executeQuery();
            while (resultSet.next()) {
                sequenceIds.add(resultSet.getLong(1));
            }
            resultSet.close();
        }
        stat.close();
        connection.close();
        return sequenceIds;
    }
    
    private List<Integer> getNeededTaxonIds() {
        List<Integer> rtns = new ArrayList<Integer>();
        rtns.add(9606); // homo sapiens
////        rtns.add(4932); // S. cerevisiae
        rtns.add(559292); // S.cerevisae S288C: Have to use this yeast. Don't use 4932!!! Otherwise, the hit is very low
        rtns.add(6239); // C. elegans
        rtns.add(7227); // D. melanogaster
        rtns.add(10090); // Mus musculus
        
        // for Zebrafish check
//        rtns.add(7955);
        
        return rtns;
    }
    
}
