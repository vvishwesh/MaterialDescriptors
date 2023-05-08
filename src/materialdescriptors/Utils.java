package materialdescriptors;

import java.io.IOException;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;


/**
 *
 * @author Vishwesh Venkatraman
 */

public class Utils
{

//------------------------------------------------------------------------------

    public static boolean isNumber(char ch)
    {
        if (Character.isDigit(ch))
                return true;
        int n = Character.compare(ch, '.');
        return n == 0;
    }

//------------------------------------------------------------------------------

    public static void pause()
    {
        System.err.println("Press a key to continue");
        try
        {
            int inchar = System.in.read();
        }
        catch (IOException ioe)
        {
            System.err.println("Error reading from user");
        }
    }

//------------------------------------------------------------------------------

    private static void normalize(ArrayList<Double> values)
    {
    	double sum = 0.;
    	int size = values.size();
        for (int i = 0; i < size; i++)
            sum += values.get(i);

        for (int i = 0; i < size; i++)
            values.set(i, values.get(i)/sum);
    }

//------------------------------------------------------------------------------

    private static double round(double value, int places)
    {
    	if (places < 0)
            throw new IllegalArgumentException();
    	BigDecimal bd = new BigDecimal(value);
    	bd = bd.setScale(places, RoundingMode.HALF_UP);
    	return bd.doubleValue();
    }
    
//------------------------------------------------------------------------------

    public static double[] convertToDoubleArray(Object[] vals)
    {
        double[] res = new double[vals.length];
        for (int i=0; i<vals.length; i++)
        {
            res[i] = Double.parseDouble((String) vals[i]);
        }
        return res;
    }

//------------------------------------------------------------------------------

    public static HashMap<String,Double> parseElementalComposition(String formula) throws Exception
    {
        Map<String,Double> element_counts = new HashMap<>();
        int len = formula.length();
        for(int i = 0; i < len;)
        {
            boolean is_group = false;
            if(formula.charAt(i) == '(')
            {
                i++;
                is_group = true;
            }
            double repeat_count = 1;
            Map<String,Double> atoms_in_group = new HashMap<>();
            do
            {
                int start = i;
                int restore_i = 0;
                String element = null;
                String restore_element = null;
                while(i < len && Character.isLetter(formula.charAt(i)))
                {
                    i++;
                    element = formula.substring(start, i);

                    String element_from_map = Constants.SYMBOLSMAP.get(element);
                    if(element_from_map == null)
                        element_from_map = Constants.SYMBOLSMAP.get(element.toUpperCase());
                    if(element_from_map != null)
                    {
                        restore_i = i;
                        restore_element = element_from_map;
                    }
                }
                if(restore_element != null)
                {
                    i = restore_i;
                    element = restore_element;
                }

                if(element == null || "".equals(element))
                {
                    throw new Exception("Parse error: could not detect an element where one was expected in formula string.");
                    //System.err.println("Remaining formula to parse: " + formula.substring(i));
                    //System.exit(-1);
                }

                start = i;
                while(i < len && isNumber(formula.charAt(i)))
                    i++;

                double count;
                try
                {
                    count = Double.parseDouble(formula.substring(start, i));
                }
                catch(NumberFormatException e)
                {
                    count = 1;
                }

                atoms_in_group.put(element, count);
                if(i < len && formula.charAt(i) == ')')
                {
                    if(!is_group)
                        System.out.println("Parse error: unmatched parenthesis detected...");

                    i++;
                    is_group = false;

                    start = i;
                    while(i < len && isNumber(formula.charAt(i)))
                        i++;

                    try
                    {
                        repeat_count = Double.parseDouble(formula.substring(start, i));
                    }
                    catch(NumberFormatException e)
                    {
                        repeat_count = 1;
                    }
                }
            } while(is_group == true);

            for(String atom_type : atoms_in_group.keySet())
            {
                double current_value = 0;
                if(element_counts.containsKey(atom_type))
                    current_value = element_counts.get(atom_type);
                element_counts.put(atom_type, current_value + atoms_in_group.get(atom_type) * repeat_count);
            }
        }

        HashMap<String, Double> elMap = new HashMap<>();
        ArrayList<Double> vals = new ArrayList<>(element_counts.values());
        ArrayList<String> keyList = new ArrayList<>(element_counts.keySet());
        normalize(vals);
        for (int i=0; i<keyList.size(); i++)
            elMap.put(Constants.ELEMENTMAP.get(keyList.get(i)), round(vals.get(i), 3));
        vals.clear();
        element_counts.clear();
        keyList.clear();

        return elMap;
    }
    
//------------------------------------------------------------------------------    

    public static int getIndexOfLargest(double[] array)
    {
        if ( array == null || array.length == 0 ) 
            return -1; // null or empty

        int largest = 0;
        for ( int i = 1; i < array.length; i++ )
            if (array[i] > array[largest]) largest = i;
        return largest; // position of the first largest found
    }

//------------------------------------------------------------------------------

    public static double getValue(String elem, String type)
    {
        double val = Double.NaN;

        OUTER:
        for (Element el : Constants.ELEMENT_PROPS)
        {
            if (el.getSymbol().equals(elem))
            {
                switch (type)
                {
                    case "EN_PAULING":
                        val = el.getEN_Pauling();
                        break OUTER;
                    case "EN_MARTYNOV":
                        val = el.getEN_MartynovBatsanov();
                        break OUTER;
                    case "EN_GHOSH":
                        val = el.getEN_Ghosh();
                        break OUTER;
                    case "EN_MULLIKEN":
                        val = el.getEN_Mulliken();
                        break OUTER;
                    case "EN_NAGLE":
                        val = el.getEN_Nagle();
                        break OUTER;    
                    case "EN_GORDY":
                        val = el.getEN_Gordy();
                        break OUTER;
                    case "EN_RAHM":
                        val = el.getEN_Rahm();
                        break OUTER;
                    case "EN_ALLEN":
                        val = el.getEN_Allen();
                        break OUTER;
                    case "EN_ALLRED":
                        val = el.getEN_AllredRochow();
                        break OUTER;
                    case "EN_COTRELL":
                        val = el.getEN_CotrellSutton();
                        break OUTER;
                    case "EN_BOEYENS":
                        val = el.getEN_Boeyens();
                        break OUTER;
                    case "ELECTROPHILICITY":
                        val = el.getElectrophilicity();
                        break OUTER;
                    case "BP":
                        val = el.getBP();
                        break OUTER;
                    case "THCOND":
                        val = el.getThermalConductivity();
                        break OUTER;
                    case "COHESIVEENG":
                        val = el.getCohesiveEnergy();
                        break OUTER;
                    case "HTFUSION":
                        val = el.getHeatOfFusion();
                        break OUTER;
                    case "HTVAP":
                        val = el.getHeatOfVaporization();
                        break OUTER;
                    case "DIPOLE":
                        val = el.getDipolePolarizability();
                        break OUTER;
                    case "ATMNUM":
                        val = el.getAtomicNumber();
                        break OUTER;
                    case "ATMWT":
                        val = el.getAtomicWeight();
                        break OUTER;
                    case "DENSITY":
                        val = el.getDensity();
                        break OUTER;
                    case "VEC":
                        val = el.getVEC();
                        break OUTER;
                    case "LQUANT":
                        val = el.getLQuantumNumber();
                        break OUTER;
                    case "ATMRAD":
                        val = el.getAtomicRadius();
                        break OUTER;
                    case "MP":
                        val = el.getMP();
                        break OUTER;
                    case "ZUNGER":
                        val = el.getRadiusZunger();
                        break OUTER;
                    case "SORB":
                        val = el.getsOrbitalElectronNumber();
                        break OUTER;
                    case "PORB":
                        val = el.getpOrbitalElectronNumber();
                        break OUTER;
                    case "DORB":
                        val = el.getdOrbitalElectronNumber();
                        break OUTER;
                    case "FORB":
                        val = el.getfOrbitalElectronNumber();
                        break OUTER;
                    case "PERIOD":
                        val = el.getPeriod();
                        break OUTER;
                    case "GROUP":
                        val = el.getGroup();
                        break OUTER;
                    case "ELAFFINITY":
                        val = el.getElectronAffinity();
                        break OUTER;
                    case "BULKMOD":
                        val = el.getBulkModulus();
                        break OUTER;
                    case "ATMVOL":
                        val = el.getAtomicVolume();
                        break OUTER;
                    case "SPECIFIC_HEAT":
                        val = el.getSpecificHeat();
                        break OUTER;
                    case "TOTAL_VALENCE_ELECTRONS":
                        val = el.getTotal_valence_electrons();
                        break OUTER;
                    case "WORK_FUNCTION":
                        val = el.getWorkFunction();
                        break OUTER;
                    case "COVALENT_RADIUS":
                        val = el.getCovalent_radius();
                        break OUTER;
                    case "CRYSTAL_RADIUS":
                        val = el.getCrystalRadius();
                        break OUTER;
                    case "IONIC_RADIUS":
                        val = el.getIonicRadius();
                        break OUTER;
                    case "ATOMENTHALPY":
                        val = el.getEnthalpyAtomization();
                        break OUTER;
                    case "ELECTRON_DENSITY":
                        val = el.getElectronDensity();
                        break OUTER;
                    case "IONIZATION_ENERGY":
                        val = el.getIonizationEnergy();
                        break OUTER;
                    case "MAGNETIC_MOMENT":
                        val = el.getMagneticMoment();
                        break OUTER;
                    case "NUM_S_UNFILLED_VALENCE_ELECTRONS":
                        val = el.getNum_s_unfilled_valence_electrons();
                        break OUTER;
                    case "NUM_F_UNFILLED_VALENCE_ELECTRONS":
                        val = el.getNum_f_unfilled_valence_electrons();
                        break OUTER;
                    case "NUM_P_UNFILLED_VALENCE_ELECTRONS":
                        val = el.getNum_p_unfilled_valence_electrons();
                        break OUTER;
                    case "NUM_D_UNFILLED_VALENCE_ELECTRONS":
                        val = el.getNum_d_unfilled_valence_electrons();
                        break OUTER;
                    case "SORBITAL_RAD":
                        val = el.getRadius_s_orbital();
                        break OUTER;
                    case "PORBITAL_RAD":
                        val = el.getRadius_p_orbital();
                        break OUTER;
                    case "DORBITAL_RAD":
                        val = el.getRadius_d_orbital();
                        break OUTER;
                    case "TOTAL_UNFILLED_ELECTRONS":
                        val = el.getTotal_unfilled_valence_electrons();
                        break OUTER;
                    case "MENDEL":
                        val = el.getMendeleevNumber();
                        break OUTER;
                    case "SGN":
                        val = el.getSpaceGroupNumber();
                        break OUTER;
                    case "ZCRITICAL":
                        val = el.getCriticalNuclearCharge();
                        break OUTER;
                    case "EN_TANDON":
                        val = el.getEN_Tandon();
                        break OUTER;
                    case "RADIUS_RAHM":
                        val = el.getAtomicRadiusRahm();
                        break OUTER;
                    case "LATTICE_CONSTANT":
                        val = el.getLatticeConstant();
                        break OUTER;
                    default:
                        break;
                }
            }
        }

        return val;
    }

//------------------------------------------------------------------------------

}
