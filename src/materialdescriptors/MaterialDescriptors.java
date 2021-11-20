package materialdescriptors;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Vishwesh Venkatraman
 */
public class MaterialDescriptors 
{
    
//------------------------------------------------------------------------------    
    
    private static boolean convertToBoolean(String value) 
    {
        boolean returnValue = false;
        if ("1".equalsIgnoreCase(value) || "yes".equalsIgnoreCase(value) || 
            "true".equalsIgnoreCase(value))
            returnValue = true;
        return returnValue;
    }
    
//------------------------------------------------------------------------------    

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) 
    {
        if(args.length < 3)
        {
            System.err.println("java -jar MaterialDescriptors.jar inputFile "
                + "header (1/0/yes/no/true/false) outputFile");
            System.err.println("inputfile: a text file containing the chemical compositions.");
            System.err.println("header: whether the input file contains a header line");
            System.err.println("outputfile: name of the file to which the descriptors will be written.");
            
            System.exit(-1);
        }
        
        String inputFile = args[0];
        boolean hasHeader = convertToBoolean(args[1]);
        String outFile = args[2];
        
        
        try
        {
            ArrayList<String> formulae = IO.readFile(inputFile, hasHeader);
            File fileToSave = new File(outFile);
                
            for (int i=0; i<formulae.size(); i++)
            {
                String compnd = formulae.get(i);
                HashMap<String,Double> elMap = Utils.parseElementalComposition(compnd);
                HashMap<String, Double> desc = new HashMap<>();
                
                DescriptorCalculator.computeDescriptors(elMap, desc);
                Descriptor cdesc = new Descriptor(compnd, desc);
                if (i == 0)
                    IO.writeDescriptor(fileToSave.getAbsolutePath(), cdesc, false);
                else
                    IO.writeDescriptor(fileToSave.getAbsolutePath(), cdesc, true);
            }
        }
        catch (Exception ex)
        {
            Logger.getLogger(MaterialDescriptors.class.getName()).log(Level.SEVERE, null, ex);
            System.exit(-1);
        }
        
        System.exit(0);
    }
    
//------------------------------------------------------------------------------    
}
