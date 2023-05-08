package materialdescriptors;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.DoubleSummaryStatistics;
import java.util.HashMap;
import java.util.Map;

// https://github.com/ademl/predict_Etot_dHf/blob/master/elemental_properties.py

/**
 *
 * @author Vishwesh Venkatraman
 */
public class DescriptorCalculator
{

//------------------------------------------------------------------------------

    /**
     * Describe the atomic packing effect in multicomponent alloys. A normalized
     * parameter of the geometric packing state should be a good candidate to
     * reveal the atomic packing instability. See doi.org/10.1016/j.scriptamat.2014.09.010
     * doi.org/10.1063/1.3587228
     * @param elemComp
     * @param descVal
     */

    public static void getAtomicSizeDifference(HashMap<String, Double> elemComp,
        HashMap<String, Double> descVal)
    {
        double rad = 0, deltarad = 0, sumrad = 0;
        double[] radvals = new double[elemComp.size()];
        double[] compvals = new double[elemComp.size()];
        int i = 0, j;
        for (Map.Entry<String, Double> entry : elemComp.entrySet())
        {
            compvals[i] = entry.getValue();
            radvals[i] = Utils.getValue(entry.getKey(), "ATMRAD");
            sumrad += compvals[i] * radvals[i];
            i++;
        }

        for (i=0; i<radvals.length; i++)
        {
            rad = radvals[i]/sumrad;
            deltarad += compvals[i] * Math.pow((1 - rad), 2.0);
            i++;
        }

        descVal.put("deltaAtomRadius", Math.sqrt(deltarad)*100);

        deltarad = 0;
        for (i=0; i<radvals.length-1; i++)
        {
            for (j=i+1; j<radvals.length; j++)
            {
                rad = (radvals[i] + radvals[j])/sumrad;
                deltarad += compvals[i] * compvals[j] * Math.pow(1 - rad, 2.0);
            }
        }

        descVal.put("varianceAtomRadius", Math.sqrt(deltarad)*100);

        // radii of the smallest atom
        double minrad = Arrays.stream(radvals).min().getAsDouble();
        // radii of the largest atom
        double maxrad = Arrays.stream(radvals).max().getAsDouble();

        double f = Math.pow(minrad + sumrad, 2);
        double omega_S = 1 - Math.sqrt((f - sumrad*sumrad)/f);
        f = Math.pow(maxrad + sumrad, 2);
        double omega_L = 1 - Math.sqrt((f - sumrad*sumrad)/f);

        descVal.put("atomicPackingMisfit", omega_S/omega_L);
    }

//------------------------------------------------------------------------------

    /**
     * Statistics of the numbers and fractions of elements in an alloy, measured
     * regardless of their type. Taken from Ward et al., Acta Materialia 159 (2018) 102-111
     * @param elemComp
     * @param descVal
     */
    private static void getStochiometricAttributes(HashMap<String, Double> elemComp,
        HashMap<String, Double> descVal)
    {
        double norm_2 = 0, norm_3 = 0, norm_5 = 0, norm_7 = 0, norm_10 = 0;
        double comp;
        for (Map.Entry<String, Double> entry : elemComp.entrySet())
        {
            comp = entry.getValue();
            norm_2 += Math.pow(comp, 2);
            norm_3 += Math.pow(comp, 3);
            norm_5 += Math.pow(comp, 5);
            norm_7 += Math.pow(comp, 7);
            norm_10 += Math.pow(comp, 10);
        }

        norm_2 = Math.pow(norm_2, 0.5);
        norm_3 = Math.pow(norm_2, 1.0/3.0);
        norm_5 = Math.pow(norm_2, 0.2);
        norm_7 = Math.pow(norm_2, 1.0/7.0);
        norm_10 = Math.pow(norm_2, 0.1);

        descVal.put("Comp_L2Norm", norm_2); // 𝐿2-norm of the element fractions
        descVal.put("Comp_L3Norm", norm_3); // 𝐿3-norm of the element fractions
        descVal.put("Comp_L5Norm", norm_5); // 𝐿5-norm of the element fractions
        descVal.put("Comp_L7Norm", norm_7); // 𝐿7-norm of the element fractions
        descVal.put("Comp_L10Norm", norm_10); // 𝐿10-norm of the element fractions
    }

//------------------------------------------------------------------------------

    /**
     * See Gao et al, dx.doi.org/10.1007/s11661-015-3105-z
     * @param elemComp
     * @param descVal
     */

    private static void getDensityMix(HashMap<String, Double> elemComp,
        HashMap<String, Double> descVal)
    {
        double compval, sumdwt = 0, sumwt = 0;
        double dval, wval;

        for (Map.Entry<String, Double> entry : elemComp.entrySet())
        {
            compval = entry.getValue();
            dval = Utils.getValue(entry.getKey(), "DENSITY");
            wval = Utils.getValue(entry.getKey(), "ATMWT");
            sumwt += compval * wval;
            sumdwt += (compval * wval)/dval;
        }
        descVal.put("alloyDensity", sumwt/sumdwt);
    }

//------------------------------------------------------------------------------

    /**
     *
     * @param elemComp
     * @param descVal
     */
    public static void getVEC(HashMap<String, Double> elemComp,
        HashMap<String, Double> descVal)
    {
        double compval, vec = 0, deltavec = 0, sumvec = 0;
        double[] vecvals = new double[elemComp.size()];
        int i = 0;
        for (Map.Entry<String, Double> entry : elemComp.entrySet())
        {
            compval = entry.getValue();
            vecvals[i] = Utils.getValue(entry.getKey(), "VEC");
            sumvec += compval * vecvals[i];
            i++;
        }
        descVal.put("frWtVEC", Math.sqrt(sumvec));


        i=0;
        for (Map.Entry<String, Double> entry : elemComp.entrySet())
        {
            compval = entry.getValue();
            vec = vecvals[i]/sumvec;
            deltavec += compval * (Math.pow((1 - vec),2.0));
            i++;
        }
        descVal.put("deltaVEC", Math.sqrt(deltavec));
    }

//------------------------------------------------------------------------------

    /**
     *
     * @param elemComp
     * @param descVal
     */
    private static void getEnthalpyOfMixing(HashMap<String, Double> elemComp,
        HashMap<String, Double> descVal)
    {
        ArrayList<Double> vals = new ArrayList<>(elemComp.values());
        ArrayList<String> keyList = new ArrayList<>(elemComp.keySet());
        String elA, elB;
        double entVal, sum = 0, conc_i, conc_j;
        for (int i=0; i<(keyList.size()-1); i++)
        {
            elA = keyList.get(i);
            conc_i = vals.get(i);
            for (int j=(i+1); j<keyList.size(); j++)
            {
                elB = keyList.get(j);
                conc_j = vals.get(j);
                entVal = Miedema.getValue(elA, elB);
                sum += entVal * conc_i * conc_j;
            }
        }
        sum *= 4;
        descVal.put("deltaHMIX", sum); // KiloJoules/mol
    }

//------------------------------------------------------------------------------

    /**
     *
     * @param elemComp
     * @param descVal
     */
    private static void getEffectiveMeltingTemperature(HashMap<String, Double> elemComp,
        HashMap<String, Double> descVal)
    {
        double mp = 0, sum = 0;

        for (Map.Entry<String, Double> entry : elemComp.entrySet())
        {
            mp = Utils.getValue(entry.getKey(), "MP") + Constants.KELVIN;
            sum += entry.getValue() * mp;
        }
        descVal.put("effMP", sum);
    }

//------------------------------------------------------------------------------


    /**
     * the average volume of atom of the alloy - 10.1016/j.jallcom.2016.09.189
     * @param elemComp
     * @param descVal
     */
    private static void getAverageVolumeAA(HashMap<String, Double> elemComp,
        HashMap<String, Double> descVal)
    {
        double bmod = 0, sum = 0, avol = 0, sum2 = 0;

        for (Map.Entry<String, Double> entry : elemComp.entrySet())
        {
            bmod = Utils.getValue(entry.getKey(), "BULKMOD");
            avol = Utils.getValue(entry.getKey(), "ATMVOL");
            sum += entry.getValue() * bmod * avol;
            sum2 += entry.getValue() * bmod;
        }
        // the average volume of atom of the alloy
        descVal.put("avgVolAlloy", sum/sum2);
    }

//------------------------------------------------------------------------------

    /**
     *
     * @param elemComp
     * @param descVal
     */

    public static void getThermodynamicParameter(HashMap<String, Double> elemComp,
        HashMap<String, Double> descVal)
    {
        double T_m = descVal.get("EFFECTIVEMP");
        double entropy = descVal.get("CONFIGENTROPY");
        double hmix = descVal.get("DELTAHMIX");
        //System.err.println(T_m + " " + entropy + " " + hmix);

        double val = (T_m * entropy)/Math.abs(hmix);
        descVal.put("SSFORMATION", val);
    }

//------------------------------------------------------------------------------

    /**
     *
     * @param elemComp
     * @param descVal
     */
    public static void getConfigurationalEntropy(HashMap<String, Double> elemComp,
        HashMap<String, Double> descVal)
    {
        double compval = 0, sum = 0;
        for (Map.Entry<String, Double> entry : elemComp.entrySet())
        {
            compval = entry.getValue();
            sum += compval * Math.log(compval);
        }
        sum *= (-1) * Constants.R;
        descVal.put("configEntropy", sum/1000); // KiloJoules/K/mol
    }


//------------------------------------------------------------------------------

    /**
     * Calculate "ionicity" of the material based on different electronegativity
     * definitions.
     * <ol>
     * <li>Maximum ionic character: Max I(x,y) between any two constituents
     * <li>Mean ionic character: sum(x<sub>i</sub> * x<sub>j</sub> * I(i,j)) where
     * x<sub>i</sub> is the fraction of element i.
     * </ol>
     * @param elemComp
     * @param descVal
     * @param enegtype
     * @param prefix
     */

    public static void getIonicCharacter(HashMap<String, Double> elemComp,
        HashMap<String, Double> descVal, String enegtype, String prefix)
    {
        double[] pvals = new double[elemComp.size()];
        double[] comp = new double[elemComp.size()];

        int i = 0;
        for (Map.Entry<String, Double> entry : elemComp.entrySet())
        {
            comp[i] = entry.getValue();
            switch (enegtype)
            {
                // electronegativity
                case "EN_PAULING":
                    pvals[i] = Utils.getValue(entry.getKey(), "EN_PAULING");
                    break;
                case "EN_MARTYNOV":
                    pvals[i] = Utils.getValue(entry.getKey(), "EN_MARTYNOV");
                    break;
                case "EN_GHOSH":
                    pvals[i] = Utils.getValue(entry.getKey(), "EN_GHOSH");
                    break;
                case "EN_MULLIKEN":
                    pvals[i] = Utils.getValue(entry.getKey(), "EN_MULLIKEN");
                    break;
                case "EN_NAGLE":
                    pvals[i] = Utils.getValue(entry.getKey(), "EN_NAGLE");
                    break;
                case "EN_GORDY":
                    pvals[i] = Utils.getValue(entry.getKey(), "EN_GORDY");
                    break;
                case "EN_RAHM":
                    pvals[i] = Utils.getValue(entry.getKey(), "EN_RAHM");
                    break;
                case "EN_ALLEN":
                    pvals[i] = Utils.getValue(entry.getKey(), "EN_ALLEN");
                    break;
                case "EN_ALLRED":
                    pvals[i] = Utils.getValue(entry.getKey(), "EN_ALLRED");
                    break;
                case "EN_BOEYENS":
                    pvals[i] = Utils.getValue(entry.getKey(), "EN_BOEYENS");
                    break;
                case "ELECTROPHILICITY":
                    pvals[i] = Utils.getValue(entry.getKey(), "ELECTROPHILICITY");
                    break;
                default:
                    break;
            }
            i++;
        }


        double eneg_diff = 0, mean_ionic_char = 0, maxval = Double.MIN_VALUE;
        for (i=0; i<(pvals.length-1); i++)
        {
            for (int j=i+1; j<pvals.length; j++)
            {
                eneg_diff = 1 - Math.exp(-0.25 * Math.pow(pvals[i] - pvals[j], 2.0));
                if (eneg_diff > maxval)
                    maxval = eneg_diff;
                mean_ionic_char += comp[i] * comp[j] * eneg_diff;
            }
        }

        descVal.put("meanic_" + prefix, mean_ionic_char);
        descVal.put("maxic_" + prefix, maxval);
    }

//------------------------------------------------------------------------------

    /**
     * Composition-weighted average number of {s,p,d,f} valence electrons of each
     * constituent element divided by average number of total valence electrons
     * @param elemComp
     * @param descVal
     * @param proptype
     * @param prefix
     */

    private static void getValenceOrbitalOccupationAttributes(HashMap<String, Double> elemComp,
        HashMap<String, Double> descVal, String proptype, String prefix)
    {
        double[] pvals = new double[elemComp.size()];
        double[] tvals = new double[elemComp.size()];
        double[] comp = new double[elemComp.size()];

        int i = 0;
        for (Map.Entry<String, Double> entry : elemComp.entrySet())
        {
            comp[i] = entry.getValue();

            tvals[i] = Utils.getValue(entry.getKey(), "TOTAL_VALENCE_ELECTRONS");

            switch (proptype)
            {
                case "NUM_S_UNFILLED_VALENCE_ELECTRONS":
                    pvals[i] = Utils.getValue(entry.getKey(), "NUM_S_UNFILLED_VALENCE_ELECTRONS");
                    break;
                case "NUM_F_UNFILLED_VALENCE_ELECTRONS":
                    pvals[i] = Utils.getValue(entry.getKey(), "NUM_F_UNFILLED_VALENCE_ELECTRONS");
                    break;
                case "NUM_P_UNFILLED_VALENCE_ELECTRONS":
                    pvals[i] = Utils.getValue(entry.getKey(), "NUM_P_UNFILLED_VALENCE_ELECTRONS");
                    break;
                case "NUM_D_UNFILLED_VALENCE_ELECTRONS":
                    pvals[i] = Utils.getValue(entry.getKey(), "NUM_D_UNFILLED_VALENCE_ELECTRONS");
                    break;
                case "SORB":
                    pvals[i] = Utils.getValue(entry.getKey(), "SORB");
                    break;
                case "PORB":
                    pvals[i] = Utils.getValue(entry.getKey(), "PORB");
                    break;
                case "DORB":
                    pvals[i] = Utils.getValue(entry.getKey(), "DORB");
                    break;
                case "FORB":
                    pvals[i] = Utils.getValue(entry.getKey(), "FORB");
                    break;
                default:
                    break;
            }
            i++;
        }

        double numerator = 0, denominator = 0;
        for (i=0; i<pvals.length; i++)
        {
            numerator += pvals[i]*comp[i];
            denominator += tvals[i]*comp[i];
        }

        double frac = numerator/denominator;

        descVal.put("frelec_" + prefix, frac);
    }

//------------------------------------------------------------------------------

    /**
     * See 10.1038/srep34256
     * @param pvals
     * @param comp
     */
    private static double[] calculateGeneralizedMeans(double[] pvals, double[] comp)
    {
        int n = pvals.length;
        double mu_geometric = 0;
        double mu_quartic_harmonic = 0, mu_cubic_harmonic = 0, mu_quadratic_harmonic = 0,
            mu_harmonic = 0;
        double mu_arithmetic = 0, mu_euclidean = 0, mu_cubic = 0, mu_quartic = 0;

        double beta = 0, sumwt = 0;
        for (int i=0; i<n; i++)
        {
            sumwt += comp[i]*comp[i];
        }
        beta = 1.0/(1.0 - sumwt);

        for (int i=0; i<n; i++)
        {
            mu_geometric += comp[i] * Math.log(pvals[i]);
            mu_quartic_harmonic += comp[i] * Math.pow(pvals[i], -4);
            mu_cubic_harmonic += comp[i] * Math.pow(pvals[i], -3);
            mu_quadratic_harmonic += comp[i] * Math.pow(pvals[i], -2);
            mu_harmonic += comp[i] * Math.pow(pvals[i], -1.0);
            mu_arithmetic += comp[i] * Math.pow(pvals[i], 1.0);
            mu_euclidean += comp[i] * Math.pow(pvals[i], 2.0);
            mu_cubic += comp[i] * Math.pow(pvals[i], 3.0);
            mu_quartic += comp[i] * Math.pow(pvals[i], 4.0);
        }


        mu_quartic_harmonic = Math.pow(mu_quartic_harmonic, -0.25);
        mu_cubic_harmonic = Math.pow(mu_quartic_harmonic, -0.33);
        mu_quadratic_harmonic = Math.pow(mu_quadratic_harmonic, -0.5);
        mu_harmonic = Math.pow(mu_quartic_harmonic, -1);
        mu_geometric = Math.exp(mu_geometric);
        mu_euclidean = Math.pow(mu_euclidean, 0.5);
        mu_cubic = Math.pow(mu_cubic, 0.33);
        mu_quartic = Math.pow(mu_quartic, 0.25);

        double sigma_0 = 0, sigma_1 = 0;
        for (int i=0; i<n; i++)
        {
            sigma_1 += comp[i] * Math.pow(pvals[i] - mu_arithmetic, 2.);
            sigma_0 += comp[i] * Math.log(Math.pow(pvals[i]/mu_geometric, 2.0));
        }

        sigma_1 = Math.pow(beta * sigma_1, 0.5);
        sigma_0 = Math.pow(Math.exp(beta * sigma_0), 0.5);

        double[] res = new double[11];
        res[0] = mu_quartic_harmonic;
        res[1] = mu_cubic_harmonic;
        res[2] = mu_quadratic_harmonic;
        res[3] = mu_harmonic;
        res[4] = mu_geometric;
        res[5] = mu_arithmetic;
        res[6] = mu_euclidean;
        res[7] = mu_cubic;
        res[8] = mu_quartic;
        res[9] = sigma_0;
        res[9] = sigma_1;

        return res;
    }


//------------------------------------------------------------------------------

    /**
     * Calculates summary statistics. Harmonic weighted mean and difference
     * See doi.org/10.1038/s41467-018-06682-4.
     * @param pvals
     * @param comp
     * @return
     */

    private static double[] calculateStatistics(double[] pvals, double[] comp)
    {
        double maxval = 0, minval = 0;
        double modeval = 0, meanval = 0;
        double maxdiff = 0, meandev = 0;

        // mode is calculated as the property value corresponding to the element
        // with the highest composition
        modeval = pvals[Utils.getIndexOfLargest(comp)];

        int n = pvals.length;
        for (int i=0; i<n; i++)
        {
            meanval += pvals[i] * comp[i];
        }

        for (int i=0; i<n; i++)
        {
            meandev += comp[i] * Math.abs(pvals[i] - meanval);
        }

        DoubleSummaryStatistics dss = Arrays.stream(pvals).summaryStatistics();
        maxval = dss.getMax();
        minval = dss.getMin();
        maxdiff = maxval - minval;

        double hmean = 0, meandiff = 0;
        for (int i=0; i<n-1; i++)
        {
            for (int j=i+1; j<n; j++)
            {
                //System.err.println(comp[i] + " " + comp[j] + " " + pvals[i] + " " + pvals[j]);
                //hmean += (comp[i] + comp[j]) * ((pvals[i] * pvals[j])/((pvals[i] + pvals[j])));
                meandiff += (comp[i] + comp[j]) * Math.abs(pvals[i] - pvals[j]);
            }
        }
        //hmean = hmean/(n-1);
        meandiff = meandiff/(n-1);


        double[] res = new double[7];
        res[0] = minval;
        res[1] = maxval;
        res[2] = maxdiff;
        res[3] = modeval;
        res[4] = meanval;
        res[5] = meandev;
        res[6] = meandiff;
        //res[7] = hmean;

        return res;
    }

//------------------------------------------------------------------------------

    /**
     * Defined as the mean, mean absolute deviation, range, minimum, maximum and
     * mode of different elemental properties. For each property, the minimum,
     * maximum, and range of the values of the properties of each element present
     * in the material is computed along with the fraction-weighted mean,
     * average deviation, and mode (i.e. the property of the most prevalent element).
     * @param elemComp
     * @param descVal
     * @param proptype
     * @param prefix
     */
    private static void getElementalPropertyStatistics(HashMap<String, Double> elemComp,
        HashMap<String, Double> descVal, String proptype, String prefix)
    {
        double[] pvals = new double[elemComp.size()];
        double[] comp = new double[elemComp.size()];

        int i = 0;
        for (Map.Entry<String, Double> entry : elemComp.entrySet())
        {
            comp[i] = entry.getValue();


            switch (proptype)
            {
                // electronegativity
                case "EN_PAULING":
                    pvals[i] = Utils.getValue(entry.getKey(), "EN_PAULING");
                    break;
                case "EN_MARTYNOV":
                    pvals[i] = Utils.getValue(entry.getKey(), "EN_MARTYNOV");
                    break;
                case "EN_GHOSH":
                    pvals[i] = Utils.getValue(entry.getKey(), "EN_GHOSH");
                    break;
                case "EN_MULLIKEN":
                    pvals[i] = Utils.getValue(entry.getKey(), "EN_MULLIKEN");
                    break;
                case "EN_NAGLE":
                    pvals[i] = Utils.getValue(entry.getKey(), "EN_NAGLE");
                    break;
                case "EN_GORDY":
                    pvals[i] = Utils.getValue(entry.getKey(), "EN_GORDY");
                    break;
                case "EN_RAHM":
                    pvals[i] = Utils.getValue(entry.getKey(), "EN_RAHM");
                    break;
                case "EN_ALLEN":
                    pvals[i] = Utils.getValue(entry.getKey(), "EN_ALLEN");
                    break;
                case "EN_ALLRED":
                    pvals[i] = Utils.getValue(entry.getKey(), "EN_ALLRED");
                    break;
                case "EN_COTRELL":
                    pvals[i] = Utils.getValue(entry.getKey(), "EN_COTRELL");
                    break;
                case "EN_BOEYENS":
                    pvals[i] = Utils.getValue(entry.getKey(), "EN_BOEYENS");
                    break;
                case "EN_TANDON":
                    pvals[i] = Utils.getValue(entry.getKey(), "EN_TANDON");
                    break;
                case "ELECTROPHILICITY":
                    pvals[i] = Utils.getValue(entry.getKey(), "ELECTROPHILICITY");
                    break;
                case "C6DP":
                    pvals[i] = Utils.getValue(entry.getKey(), "C6DP");
                    break;
                // atomic number
                case "ATMNUM":
                    pvals[i] = Utils.getValue(entry.getKey(), "ATMNUM");
                    break;
                // atomic radius
                case "ATMRAD":
                    pvals[i] = Utils.getValue(entry.getKey(), "ATMRAD");
                    break;
                case "ZUNGER":
                    pvals[i] = Utils.getValue(entry.getKey(), "ZUNGER");
                    break;
                case "COVALENT_RADIUS":
                    pvals[i] = Utils.getValue(entry.getKey(), "COVALENT_RADIUS");
                    break;
                case "COVALENT_RADIUS_PYYKKO":
                    pvals[i] = Utils.getValue(entry.getKey(), "COVALENT_RADIUS_PYYKKO");
                    break;
                case "COVALENT_RADIUS_CORDERO":
                    pvals[i] = Utils.getValue(entry.getKey(), "COVALENT_RADIUS_CORDERO");
                    break;
                case "IONIC_RADIUS":
                    pvals[i] = Utils.getValue(entry.getKey(), "IONIC_RADIUS");
                    break;
                case "CRYSTAL_RADIUS":
                    pvals[i] = Utils.getValue(entry.getKey(), "CRYSTAL_RADIUS");
                    break;
                case "RADIUS_RAHM":
                    pvals[i] = Utils.getValue(entry.getKey(), "RADIUS_RAHM");
                    break;
                // atomic weight
                case "ATMWT":
                    pvals[i] = Utils.getValue(entry.getKey(), "ATMWT");
                    break;
                case "MENDEL":
                    pvals[i] = Utils.getValue(entry.getKey(), "MENDEL");
                    break;
                case "PETTIFOR":
                    pvals[i] = Utils.getValue(entry.getKey(), "PETTIFOR");
                    break;
                case "GLAWE":
                    pvals[i] = Utils.getValue(entry.getKey(), "GLAWE");
                    break;
                case "DIPOLE":
                    pvals[i] = Utils.getValue(entry.getKey(), "DIPOLE");
                    break;
                case "POLARIZABILITY":
                    pvals[i] = Utils.getValue(entry.getKey(), "POLARIZABILITY");
                    break;
                case "COHESIVEENG":
                    pvals[i] = Utils.getValue(entry.getKey(), "COHESIVEENG");
                    break;
                case "HTVAP":
                    pvals[i] = Utils.getValue(entry.getKey(), "HTVAP");
                    break;
                case "ATOMENTHALPY":
                    pvals[i] = Utils.getValue(entry.getKey(), "ATOMENTHALPY");
                    break;
                case "HTFUSION":
                    pvals[i] = Utils.getValue(entry.getKey(), "HTFUSION");
                    break;
                case "HOF":
                    pvals[i] = Utils.getValue(entry.getKey(), "HOF");
                    break;
                case "SPECIFIC_HEAT":
                    pvals[i] = Utils.getValue(entry.getKey(), "SPECIFIC_HEAT");
                    break;
                case "ELECTRON_DENSITY":
                    pvals[i] = Utils.getValue(entry.getKey(), "ELECTRON_DENSITY");
                    break;
                case "ZCRITICAL":
                    pvals[i] = Utils.getValue(entry.getKey(), "ZCRITICAL");
                    break;
                case "ELAFFINITY":
                    pvals[i] = Utils.getValue(entry.getKey(), "ELAFFINITY");
                    break;
                case "MAGNETIC_MOMENT":
                    pvals[i] = Utils.getValue(entry.getKey(), "MAGNETIC_MOMENT");
                    break;
                case "IONIZATION_ENERGY":
                    pvals[i] = Utils.getValue(entry.getKey(), "IONIZATION_ENERGY");
                    break;
                case "SGN":
                    pvals[i] = Utils.getValue(entry.getKey(), "SGN");
                    break;
                case "WORK_FUNCTION":
                    pvals[i] = Utils.getValue(entry.getKey(), "WORK_FUNCTION");
                    break;
                case "THCOND":
                    pvals[i] = Utils.getValue(entry.getKey(), "THCOND");
                    break;
                case "MP":
                    pvals[i] = Utils.getValue(entry.getKey(), "MP");
                    break;
                case "BP":
                    pvals[i] = Utils.getValue(entry.getKey(), "BP");
                    break;
                // valence/orbital based
                case "LQUANT":
                    pvals[i] = Utils.getValue(entry.getKey(), "LQUANT");
                    break;
                case "TOTAL_VALENCE_ELECTRONS":
                    pvals[i] = Utils.getValue(entry.getKey(), "TOTAL_VALENCE_ELECTRONS");
                    break;
                case "NUM_S_UNFILLED_VALENCE_ELECTRONS":
                    pvals[i] = Utils.getValue(entry.getKey(), "NUM_S_UNFILLED_VALENCE_ELECTRONS");
                    break;
                case "NUM_F_UNFILLED_VALENCE_ELECTRONS":
                    pvals[i] = Utils.getValue(entry.getKey(), "NUM_F_UNFILLED_VALENCE_ELECTRONS");
                    break;
                case "NUM_P_UNFILLED_VALENCE_ELECTRONS":
                    pvals[i] = Utils.getValue(entry.getKey(), "NUM_P_UNFILLED_VALENCE_ELECTRONS");
                    break;
                case "NUM_D_UNFILLED_VALENCE_ELECTRONS":
                    pvals[i] = Utils.getValue(entry.getKey(), "NUM_D_UNFILLED_VALENCE_ELECTRONS");
                    break;
                case "SORBITAL_RAD":
                    pvals[i] = Utils.getValue(entry.getKey(), "SORBITAL_RAD");
                    break;
                case "PORBITAL_RAD":
                    pvals[i] = Utils.getValue(entry.getKey(), "PORBITAL_RAD");
                    break;
                case "DORBITAL_RAD":
                    pvals[i] = Utils.getValue(entry.getKey(), "DORBITAL_RAD");
                    break;
                case "SORB":
                    pvals[i] = Utils.getValue(entry.getKey(), "SORB");
                    break;
                case "PORB":
                    pvals[i] = Utils.getValue(entry.getKey(), "PORB");
                    break;
                case "DORB":
                    pvals[i] = Utils.getValue(entry.getKey(), "DORB");
                    break;
                case "FORB":
                    pvals[i] = Utils.getValue(entry.getKey(), "FORB");
                    break;
                case "TOTAL_UNFILLED_ELECTRONS":
                    pvals[i] = Utils.getValue(entry.getKey(), "TOTAL_UNFILLED_ELECTRONS");
                    break;
                case "LATTICE_CONSTANT":
                    pvals[i] = Utils.getValue(entry.getKey(), "LATTICE_CONSTANT");
                    break;
                case "GROUP":
                    pvals[i] = Utils.getValue(entry.getKey(), "GROUP");
                    break;
                case "PERIOD":
                    pvals[i] = Utils.getValue(entry.getKey(), "PERIOD");
                    break;
                default:
                    break;
            }

            i++;
        }

        double[] res = calculateStatistics(pvals, comp);
        descVal.put("min_" + prefix, res[0]);
        descVal.put("max_" + prefix, res[1]);
        descVal.put("maxdiff_" + prefix, res[2]);
        descVal.put("mode_" + prefix, res[3]);
        descVal.put("fwtmean_" + prefix, res[4]); // fraction weighted
        descVal.put("fwtmeandev_" + prefix, res[5]); // fraction weighted
        descVal.put("meandiff_" + prefix, res[6]); // stoichiometrically weighted mean difference
    }

//------------------------------------------------------------------------------

    private static void getHolderMeans(HashMap<String, Double> elemComp,
        HashMap<String, Double> descVal, String proptype, String prefix)
    {
        double[] pvals = new double[elemComp.size()];
        double[] comp = new double[elemComp.size()];

        int i = 0;
        for (Map.Entry<String, Double> entry : elemComp.entrySet())
        {
            comp[i] = entry.getValue();

            switch (proptype)
            {
                // electronegativity
                case "EN_PAULING":
                    pvals[i] = Utils.getValue(entry.getKey(), "EN_PAULING");
                    break;
                case "EN_MARTYNOV":
                    pvals[i] = Utils.getValue(entry.getKey(), "EN_MARTYNOV");
                    break;
                case "EN_GHOSH":
                    pvals[i] = Utils.getValue(entry.getKey(), "EN_GHOSH");
                    break;
                case "EN_MULLIKEN":
                    pvals[i] = Utils.getValue(entry.getKey(), "EN_MULLIKEN");
                    break;
                case "EN_NAGLE":
                    pvals[i] = Utils.getValue(entry.getKey(), "EN_NAGLE");
                    break;
                case "EN_GORDY":
                    pvals[i] = Utils.getValue(entry.getKey(), "EN_GORDY");
                    break;
                case "EN_RAHM":
                    pvals[i] = Utils.getValue(entry.getKey(), "EN_RAHM");
                    break;
                case "EN_ALLEN":
                    pvals[i] = Utils.getValue(entry.getKey(), "EN_ALLEN");
                    break;
                case "EN_ALLRED":
                    pvals[i] = Utils.getValue(entry.getKey(), "EN_ALLRED");
                    break;
                case "EN_COTRELL":
                    pvals[i] = Utils.getValue(entry.getKey(), "EN_COTRELL");
                    break;
                case "EN_BOEYENS":
                    pvals[i] = Utils.getValue(entry.getKey(), "EN_BOEYENS");
                    break;
                case "EN_TANDON":
                    pvals[i] = Utils.getValue(entry.getKey(), "EN_TANDON");
                    break;
                case "ELECTROPHILICITY":
                    pvals[i] = Utils.getValue(entry.getKey(), "ELECTROPHILICITY");
                    break;
                case "C6DP":
                    pvals[i] = Utils.getValue(entry.getKey(), "C6DP");
                    break;
                // atomic number
                case "ATMNUM":
                    pvals[i] = Utils.getValue(entry.getKey(), "ATMNUM");
                    break;
                // atomic radius
                case "ATMRAD":
                    pvals[i] = Utils.getValue(entry.getKey(), "ATMRAD");
                    break;
                case "ZUNGER":
                    pvals[i] = Utils.getValue(entry.getKey(), "ZUNGER");
                    break;
                case "COVALENT_RADIUS":
                    pvals[i] = Utils.getValue(entry.getKey(), "COVALENT_RADIUS");
                    break;
                case "COVALENT_RADIUS_PYYKKO":
                    pvals[i] = Utils.getValue(entry.getKey(), "COVALENT_RADIUS_PYYKKO");
                    break;
                case "COVALENT_RADIUS_CORDERO":
                    pvals[i] = Utils.getValue(entry.getKey(), "COVALENT_RADIUS_CORDERO");
                    break;
                case "IONIC_RADIUS":
                    pvals[i] = Utils.getValue(entry.getKey(), "IONIC_RADIUS");
                    break;
                case "CRYSTAL_RADIUS":
                    pvals[i] = Utils.getValue(entry.getKey(), "CRYSTAL_RADIUS");
                    break;
                case "RADIUS_RAHM":
                    pvals[i] = Utils.getValue(entry.getKey(), "RADIUS_RAHM");
                    break;
                // atomic weight
                case "ATMWT":
                    pvals[i] = Utils.getValue(entry.getKey(), "ATMWT");
                    break;
                case "MENDEL":
                    pvals[i] = Utils.getValue(entry.getKey(), "MENDEL");
                    break;
                case "PETTIFOR":
                    pvals[i] = Utils.getValue(entry.getKey(), "PETTIFOR");
                    break;
                case "GLAWE":
                    pvals[i] = Utils.getValue(entry.getKey(), "GLAWE");
                    break;
                case "DIPOLE":
                    pvals[i] = Utils.getValue(entry.getKey(), "DIPOLE");
                    break;
                case "POLARIZABILITY":
                    pvals[i] = Utils.getValue(entry.getKey(), "POLARIZABILITY");
                    break;
                case "COHESIVEENG":
                    pvals[i] = Utils.getValue(entry.getKey(), "COHESIVEENG");
                    break;
                case "HTVAP":
                    pvals[i] = Utils.getValue(entry.getKey(), "HTVAP");
                    break;
                case "ATOMENTHALPY":
                    pvals[i] = Utils.getValue(entry.getKey(), "ATOMENTHALPY");
                    break;
                case "HTFUSION":
                    pvals[i] = Utils.getValue(entry.getKey(), "HTFUSION");
                    break;
                case "HOF":
                    pvals[i] = Utils.getValue(entry.getKey(), "HOF");
                    break;
                case "SPECIFIC_HEAT":
                    pvals[i] = Utils.getValue(entry.getKey(), "SPECIFIC_HEAT");
                    break;
                case "ELECTRON_DENSITY":
                    pvals[i] = Utils.getValue(entry.getKey(), "ELECTRON_DENSITY");
                    break;
                case "ZCRITICAL":
                    pvals[i] = Utils.getValue(entry.getKey(), "ZCRITICAL");
                    break;
                case "ELAFFINITY":
                    pvals[i] = Utils.getValue(entry.getKey(), "ELAFFINITY");
                    break;
                case "MAGNETIC_MOMENT":
                    pvals[i] = Utils.getValue(entry.getKey(), "MAGNETIC_MOMENT");
                    break;
                case "IONIZATION_ENERGY":
                    pvals[i] = Utils.getValue(entry.getKey(), "IONIZATION_ENERGY");
                    break;
                case "SGN":
                    pvals[i] = Utils.getValue(entry.getKey(), "SGN");
                    break;
                case "WORK_FUNCTION":
                    pvals[i] = Utils.getValue(entry.getKey(), "WORK_FUNCTION");
                    break;
                case "THCOND":
                    pvals[i] = Utils.getValue(entry.getKey(), "THCOND");
                    break;
                case "MP":
                    pvals[i] = Utils.getValue(entry.getKey(), "MP");
                    break;
                case "BP":
                    pvals[i] = Utils.getValue(entry.getKey(), "BP");
                    break;
                // valence/orbital based
                case "LQUANT":
                    pvals[i] = Utils.getValue(entry.getKey(), "LQUANT");
                    break;
                case "TOTAL_VALENCE_ELECTRONS":
                    pvals[i] = Utils.getValue(entry.getKey(), "TOTAL_VALENCE_ELECTRONS");
                    break;
                case "NUM_S_UNFILLED_VALENCE_ELECTRONS":
                    pvals[i] = Utils.getValue(entry.getKey(), "NUM_S_UNFILLED_VALENCE_ELECTRONS");
                    break;
                case "NUM_F_UNFILLED_VALENCE_ELECTRONS":
                    pvals[i] = Utils.getValue(entry.getKey(), "NUM_F_UNFILLED_VALENCE_ELECTRONS");
                    break;
                case "NUM_P_UNFILLED_VALENCE_ELECTRONS":
                    pvals[i] = Utils.getValue(entry.getKey(), "NUM_P_UNFILLED_VALENCE_ELECTRONS");
                    break;
                case "NUM_D_UNFILLED_VALENCE_ELECTRONS":
                    pvals[i] = Utils.getValue(entry.getKey(), "NUM_D_UNFILLED_VALENCE_ELECTRONS");
                    break;
                case "SORBITAL_RAD":
                    pvals[i] = Utils.getValue(entry.getKey(), "SORBITAL_RAD");
                    break;
                case "PORBITAL_RAD":
                    pvals[i] = Utils.getValue(entry.getKey(), "PORBITAL_RAD");
                    break;
                case "DORBITAL_RAD":
                    pvals[i] = Utils.getValue(entry.getKey(), "DORBITAL_RAD");
                    break;
                case "SORB":
                    pvals[i] = Utils.getValue(entry.getKey(), "SORB");
                    break;
                case "PORB":
                    pvals[i] = Utils.getValue(entry.getKey(), "PORB");
                    break;
                case "DORB":
                    pvals[i] = Utils.getValue(entry.getKey(), "DORB");
                    break;
                case "FORB":
                    pvals[i] = Utils.getValue(entry.getKey(), "FORB");
                    break;
                case "TOTAL_UNFILLED_ELECTRONS":
                    pvals[i] = Utils.getValue(entry.getKey(), "TOTAL_UNFILLED_ELECTRONS");
                    break;
                case "LATTICE_CONSTANT":
                    pvals[i] = Utils.getValue(entry.getKey(), "LATTICE_CONSTANT");
                    break;
                case "GROUP":
                    pvals[i] = Utils.getValue(entry.getKey(), "GROUP");
                    break;
                case "PERIOD":
                    pvals[i] = Utils.getValue(entry.getKey(), "PERIOD");
                    break;
                default:
                    break;
            }

            i++;
        }

        double[] res = calculateGeneralizedMeans(pvals, comp);
        descVal.put("mu_quartic_harmonic_" + prefix, res[0]);
        descVal.put("mu_cubic_harmonic_" + prefix, res[1]);
        descVal.put("mu_quadratic_harmonic_" + prefix, res[2]);
        descVal.put("mu_harmonic_" + prefix, res[3]);
        descVal.put("mu_geometric_" + prefix, res[4]);
        descVal.put("mu_arithmetic_" + prefix, res[5]);
        descVal.put("mu_euclidean_" + prefix, res[6]);
        descVal.put("mu_cubic_" + prefix, res[7]);
        descVal.put("mu_quartic_" + prefix, res[8]);
        descVal.put("sd_geometric_" + prefix, res[9]);
        descVal.put("sd_arithmetic_" + prefix, res[10]);

    }

//------------------------------------------------------------------------------

    public static void computeStandardSet(HashMap<String, Double> elemComp,
        HashMap<String, Double> descVal)
    {
        getAverageVolumeAA(elemComp, descVal);
        getDensityMix(elemComp, descVal);
        getStochiometricAttributes(elemComp, descVal);
        getConfigurationalEntropy(elemComp, descVal);
        getEffectiveMeltingTemperature(elemComp, descVal);
        getAtomicSizeDifference(elemComp, descVal);


        getElementalPropertyStatistics(elemComp, descVal, "EN_RAHM", "eneg_rahm");
        getElementalPropertyStatistics(elemComp, descVal, "NUM_S_UNFILLED_VALENCE_ELECTRONS", "NsUnfValence");
        getElementalPropertyStatistics(elemComp, descVal, "NUM_P_UNFILLED_VALENCE_ELECTRONS", "NpUnfValence");
        getElementalPropertyStatistics(elemComp, descVal, "NUM_F_UNFILLED_VALENCE_ELECTRONS", "NfUnfValence");
        getElementalPropertyStatistics(elemComp, descVal, "NUM_D_UNFILLED_VALENCE_ELECTRONS", "NdUnfValence");
        getElementalPropertyStatistics(elemComp, descVal, "SORB", "NsValence");
        getElementalPropertyStatistics(elemComp, descVal, "PORB", "NpValence");
        getElementalPropertyStatistics(elemComp, descVal, "FORB", "NfValence");
        getElementalPropertyStatistics(elemComp, descVal, "DORB", "NdValence");


        getValenceOrbitalOccupationAttributes(elemComp, descVal, "NUM_S_UNFILLED_VALENCE_ELECTRONS", "NsUnfValence");
        getValenceOrbitalOccupationAttributes(elemComp, descVal, "NUM_P_UNFILLED_VALENCE_ELECTRONS", "NpUnfValence");
        getValenceOrbitalOccupationAttributes(elemComp, descVal, "NUM_F_UNFILLED_VALENCE_ELECTRONS", "NfUnfValence");
        getValenceOrbitalOccupationAttributes(elemComp, descVal, "NUM_D_UNFILLED_VALENCE_ELECTRONS", "NdUnfValence");
        getValenceOrbitalOccupationAttributes(elemComp, descVal, "SORB", "NsValence");
        getValenceOrbitalOccupationAttributes(elemComp, descVal, "PORB", "NpValence");
        getValenceOrbitalOccupationAttributes(elemComp, descVal, "FORB", "NfValence");
        getValenceOrbitalOccupationAttributes(elemComp, descVal, "DORB", "NdValence");

        getElementalPropertyStatistics(elemComp, descVal, "ATMNUM", "atmnum");
        getElementalPropertyStatistics(elemComp, descVal, "ATMWT", "atmwt");
        getElementalPropertyStatistics(elemComp, descVal, "RADIUS_RAHM", "rahmrad");
        getElementalPropertyStatistics(elemComp, descVal, "LQUANT", "lquant");
        getElementalPropertyStatistics(elemComp, descVal, "ZUNGER", "zungerad");
        getElementalPropertyStatistics(elemComp, descVal, "LATTICE_CONSTANT", "latconst");
        getElementalPropertyStatistics(elemComp, descVal, "MENDEL", "mendeleevnum");
        getElementalPropertyStatistics(elemComp, descVal, "DIPOLE", "dipole");
        getElementalPropertyStatistics(elemComp, descVal, "IONIZATION_ENERGY", "ioneng");
        getElementalPropertyStatistics(elemComp, descVal, "COHESIVEENG", "cohesiveng");
        getElementalPropertyStatistics(elemComp, descVal, "ELAFFINITY", "elaff");
        getElementalPropertyStatistics(elemComp, descVal, "ATOMENTHALPY", "enthalpyAtomization");
        getElementalPropertyStatistics(elemComp, descVal, "MP", "mp");
        getElementalPropertyStatistics(elemComp, descVal, "BP", "bp");
        getElementalPropertyStatistics(elemComp, descVal, "HTFUSION", "htfusion");
        getElementalPropertyStatistics(elemComp, descVal, "HTVAP", "htevap");
        getElementalPropertyStatistics(elemComp, descVal, "PERIOD", "period");
        getElementalPropertyStatistics(elemComp, descVal, "WORK_FUNCTION", "workfn");
        getElementalPropertyStatistics(elemComp, descVal, "ELECTRON_DENSITY", "eden");
        getElementalPropertyStatistics(elemComp, descVal, "SPECIFIC_HEAT", "spheat");
        getElementalPropertyStatistics(elemComp, descVal, "SGN", "sgn");
        getElementalPropertyStatistics(elemComp, descVal, "ZCRITICAL", "critnuccharge");
        getElementalPropertyStatistics(elemComp, descVal, "TOTAL_VALENCE_ELECTRONS", "totvalelec");
        getElementalPropertyStatistics(elemComp, descVal, "TOTAL_UNFILLED_ELECTRONS", "totunfilledelec");
        getElementalPropertyStatistics(elemComp, descVal, "DORBITAL_RAD", "rad_d_orb");
        getElementalPropertyStatistics(elemComp, descVal, "PORBITAL_RAD", "rad_p_orb");
        getElementalPropertyStatistics(elemComp, descVal, "SORBITAL_RAD", "rad_s_orb");
    }

//------------------------------------------------------------------------------

}
