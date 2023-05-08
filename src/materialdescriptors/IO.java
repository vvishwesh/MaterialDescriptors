package materialdescriptors;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.LineIterator;

/**
 *
 * @author Vishwesh Venkatraman
 */
public class IO
{

//------------------------------------------------------------------------------

    public static ArrayList<String> readFile(String filename, boolean header) throws IOException, Exception
    {
        ArrayList<String> fdata = new ArrayList<>();
        File file = new File(filename);
        LineIterator iterator = null;
        int ln = 0;
        try
        {
            iterator = FileUtils.lineIterator(file, "UTF-8");
            while (iterator.hasNext())
            {
                String data = iterator.nextLine();
                ln++;
                if (ln == 1 && header == true)
                {
                }
                else
                {
                    String[] arr = data.split("\\s+");
                    fdata.add(arr[0]);
                }
            }
        }
        catch (FileNotFoundException ex)
        {
            System.err.println(ex.getMessage());
            throw ex;
        }
        catch (IOException ex)
        {
            System.err.println(ex.getMessage());
            throw ex;
        }
        finally
        {
            LineIterator.closeQuietly(iterator);
        }

        return fdata;
    }

//------------------------------------------------------------------------------

    public static void writeDescriptor(String filename, Descriptor cdesc,
        boolean append) throws Exception
    {
        StringBuilder sbHeader = new StringBuilder(256);
        StringBuilder sbDesc = new StringBuilder(1024);

        
        if (!append)
        {
            ArrayList<String> keyList = new ArrayList<>(cdesc.getDescriptorValues().keySet());

            for (String str:keyList)
            {
                sbHeader.append(str).append(" ");
            }
        }

        sbDesc.append(cdesc.getCompoundName()).append(" ");


        HashMap<String, Double> map = cdesc.getDescriptorValues();


        for (Map.Entry<String, Double> entry : map.entrySet())
        {
            //System.out.println("Key = " + entry.getKey() + ", Value = " + entry.getValue());

            if (Double.isNaN(entry.getValue()) || Double.isInfinite(entry.getValue()))
            {
                sbDesc.append("NA").append(" ");
            }
            else
                sbDesc.append(String.format("%.3f", entry.getValue())).append(" ");
        }
        sbDesc.append(System.getProperty("line.separator"));

        FileWriter fw = null;
        try
        {
            fw = new FileWriter(new File(filename), append);
            if (sbHeader.toString().trim().length() > 0)
                fw.write(sbHeader.toString().trim() + System.getProperty("line.separator"));
            fw.write(sbDesc.toString().trim() + System.getProperty("line.separator"));
            fw.flush();
        }
        catch (IOException ioe)
        {
            throw ioe;
        }
        finally
        {
            sbHeader.setLength(0);
            sbDesc.setLength(0);
            try
            {
                if (fw != null)
                {
                    fw.close();
                }
            }
            catch (IOException ioe)
            {
                throw ioe;
            }
        }

    }

//------------------------------------------------------------------------------

    public static void writeDescriptors(String filename, ArrayList<Descriptor> lstDesc,
        boolean append) throws Exception
    {
        StringBuilder sbHeader = new StringBuilder(256);
        StringBuilder sbDesc = new StringBuilder(1024);


        boolean hasNA = false;

        Descriptor cdesc = lstDesc.get(0);
        ArrayList<String> keyList = new ArrayList<>(cdesc.getDescriptorValues().keySet());


        for (String str:keyList)
        {
            sbHeader.append(str).append(" ");
        }

        for (int i=0; i<lstDesc.size(); i++)
        {
            cdesc = lstDesc.get(i);
            sbDesc.append(cdesc.getCompoundName()).append(" ");
            HashMap<String, Double> map = cdesc.getDescriptorValues();

            for (Map.Entry<String, Double> entry : map.entrySet())
            {
                //System.out.println("Key = " + entry.getKey() + ", Value = " + entry.getValue());

                if (Double.isNaN(entry.getValue()) || Double.isInfinite(entry.getValue()))
                {
                    sbDesc.append("NA").append(" ");
                    hasNA = true;
                }
                else
                    sbDesc.append(String.format("%.3f", entry.getValue())).append(" ");
            }
            sbDesc.append(System.getProperty("line.separator"));
        }


        FileWriter fw = null;
        try
        {
            fw = new FileWriter(new File(filename), append);
            fw.write(sbHeader.toString().trim() + System.getProperty("line.separator"));
            fw.write(sbDesc.toString().trim() + System.getProperty("line.separator"));
            fw.flush();
        }
        catch (IOException ioe)
        {
            throw ioe;
        }
        finally
        {
            sbHeader.setLength(0);
            sbDesc.setLength(0);
            try
            {
                if (fw != null)
                {
                    fw.close();
                }
            }
            catch (IOException ioe)
            {
                throw ioe;
            }
        }

        if (hasNA)
        {
            System.err.println("#---------------------------------------#");
            System.err.println("Warning: some descriptors have NA values.");
            System.err.println("#---------------------------------------#");
        }
    }

//------------------------------------------------------------------------------

}

