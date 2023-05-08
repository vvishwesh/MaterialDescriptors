package materialdescriptors;

import java.util.ArrayList;
import java.util.HashMap;

/**
 *
 * @author Vishwesh Venkatraman
 */

public class Descriptor
{
    private String compound;
    private HashMap<String, Double> descVal;
    
//------------------------------------------------------------------------------    
    
    public Descriptor(String cmpd, HashMap<String, Double> vals)
    {
        compound = cmpd;
        descVal = vals;
    }
    
//------------------------------------------------------------------------------    

    public String getCompoundName()
    {
        return compound;
    }
    
//------------------------------------------------------------------------------    

    public void setCompoundName(String compound)
    {
        this.compound = compound;
    }
    
//------------------------------------------------------------------------------    

    public HashMap<String, Double> getDescriptorValues()
    {
        return descVal;
    }
    
//------------------------------------------------------------------------------    

    public void setDescriptorValues(HashMap<String, Double> descVal)
    {
        this.descVal = descVal;
    }
    
//------------------------------------------------------------------------------    
    
    public void cleanup()
    {
        this.descVal.clear();
    }
    
//------------------------------------------------------------------------------            
    
    public String getStringRepresentation()
    {
        StringBuilder sb = new StringBuilder(64);
        ArrayList<String> keyList = new ArrayList<>(descVal.keySet());
        sb.append("COMPOUND").append(" ");
        for (int i=0; i<keyList.size(); i++)
        {
            sb.append(keyList.get(i)).append(" ");
        }
        sb.append(System.getProperty("line.separator"));
        sb.append(compound).append(" ");
        ArrayList<Double> valList = new ArrayList<>(descVal.values());
        for (int i=0; i<valList.size(); i++)
        {
            sb.append(valList.get(i)).append(" ");
        }
        sb.append(System.getProperty("line.separator"));
        keyList.clear();
        valList.clear();
        return sb.toString().trim();
    }
    
//------------------------------------------------------------------------------                
}

