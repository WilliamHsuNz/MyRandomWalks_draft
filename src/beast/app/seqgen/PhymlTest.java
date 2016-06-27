import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;
import java.util.ArrayList;

class PhymlTest{

    public void test() {
        String phymlOutput = null;
        BufferedReader br = null;
        try {
            br = new BufferedReader(new FileReader("/Users/williamhsu/Documents/2016 Computer Science/Compsci 789 A&B/code/Phyml/output.txt"));
            StringBuilder sb = new StringBuilder();
            String line = br.readLine();
            while (line != null) {
                sb.append(line);
                sb.append(System.lineSeparator());
                line = br.readLine();
            }
            phymlOutput = sb.toString();
        } catch (IOException e) {
            System.err.println("Caught IOException: " + e.getMessage());
        } finally {
            try {
                br.close();
            } catch (IOException e) {
                System.err.println("Caught IOException2: " + e.getMessage());
            }
        }
        //System.out.println(phymlOutput);
        String[] phymlOutputArray = phymlOutput.split("\\s+");
        System.out.println(phymlOutputArray.length);
        //for(String value: phymlOutputArray){
        //	System.out.println(value);
        //}
        String tree = phymlOutputArray[52]; // get newick tree
        System.out.print(phymlOutputArray[52]);
        String[] phymlArray = new String[599];
        System.arraycopy(phymlOutputArray, 57, phymlArray, 0, 599);
        HashMap<String, double[]> ttmValues = new HashMap<String, double[]>();
        for (int i = 0; i < 599; i = i + 6) {
            String tLab = phymlArray[i].substring(0,5);
            double lat = Double.parseDouble(phymlArray[i + 2]);
            double lon = Double.parseDouble(phymlArray[i + 4].substring(0, 8));

            System.out.println(tLab +" : "+ String.valueOf(lat) + ", "+ String.valueOf(lon));
            double[]latLon = {lat, lon};
            ttmValues.put(tLab, latLon);
        }
    }


    public static void main(String [] args){
        PhymlTest pTest = new PhymlTest();
        pTest.test();
    }
}