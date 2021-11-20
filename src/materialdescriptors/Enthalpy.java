package materialdescriptors;

/**
 *
 * @author Vishwesh Venkatraman
 */
public class Enthalpy
{
    private String elemA;
    private String elemB;
    private double value;
    
    public Enthalpy(String eA, String eB, double val)
    {
        elemA = eA;
        elemB = eB;
        value = val;
    }

    public String getElemA()
    {
        return elemA;
    }

    public void setElemA(String elemA)
    {
        this.elemA = elemA;
    }

    public String getElemB()
    {
        return elemB;
    }

    public void setElemB(String elemB)
    {
        this.elemB = elemB;
    }

    public double getValue()
    {
        return value;
    }

    public void setValue(double value)
    {
        this.value = value;
    }   
}
