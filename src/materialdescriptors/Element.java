package materialdescriptors;

/**
 *
 * @author Vishwesh Venkatraman
 */
public class Element
{
    private String symbol;
    private int atomicNumber;
    private double atomicWeight;
    private double VEC; // valence electron concentration
    private double atomicRadius;
    private double vdwRadius;    
    private double MP;
    private double BP;
    private double thermalConductivity;
    private double cohesiveEnergy;
    private double ionizationEnergy;
    private double dipolePolarizability;
    private double heatOfFusion;
    private double heatOfVaporization;
    private int mendeleevNumber;
    private int LQuantumNumber;
    private double density;
    private double electronAffinity;
    private int sOrbitalElectronNumber;
    private int pOrbitalElectronNumber;
    private int dOrbitalElectronNumber;
    private int fOrbitalElectronNumber;
    private int valence;
    private double group;
    private double period;
    private int spaceGroupNumber;
    private double specificHeat;
    
    private double radiusZunger;
    
    private double atomicVolume;
    private double covalent_radius;
    private double bulkModulus;
    private double electronDensity;
    private double magneticMoment;
    private double workFunction;
    private int num_f_unfilled_valence_electrons;
    private int num_p_unfilled_valence_electrons;
    private int num_s_unfilled_valence_electrons;
    private int num_d_unfilled_valence_electrons;
    private double radius_d_orbital;
    private double radius_p_orbital;
    private double radius_s_orbital;
    private int total_unfilled_valence_electrons;
    private int total_valence_electrons;

    private double en_MartynovBatsanov;
    private double en_Pauling;
    private double en_Ghosh;
    private double en_Mulliken;
    private double en_Nagle;
    private double en_Gordy;
    private double en_Rahm;
    private double en_Allen;
    private double en_AllredRochow;
    private double en_CotrellSutton;
    private double en_Boeyens;
    private double electrophilicity;
    private double ionicRadius;
    private double crystalRadius;
    private double enthalpyAtomization;
    private double critical_nuclear_charge;
    private double en_Tandon;
    private double latticeConstant;
    private double atomicRadiusRahm;
    

    public Element(String _symbol, int _atomicNumber, double _atomicWeight, double _atomicRadius, 
        int _mendeleevNumber, double _vdwRadius, double _electronegativityPauling, int _LQuantumNumber,
        int _sOrbital, int _pOrbital, int _dOrbital, int _fOrbital,
        double _ionizationEnergy, double _electronAffinity, double _dipolePolarizability, double _heatOfFusion,
        double _heatOfVaporization, double _thermalConductivity, double _cohesiveEnergy,
        double _VEC, double _MP, double _density, double _BP, double _group, double _period, 
        int _SGN, double _electronegativityMB, double _zungerRadius, 
        double _atomicVolume, double _covalent_radius, double _bulkModulus, 
        double _electronDensity, double _magneticMoment, double _workFunction, 
        int _num_f_unfilled_valence_electrons, int _num_p_unfilled_valence_electrons, 
        int _num_s_unfilled_valence_electrons, int _num_d_unfilled_valence_electrons,
        double _radius_d_orbital, double _radius_p_orbital, double _radius_s_orbital,         
        int _total_unfilled_valence_electrons, int _total_valence_electrons, 
        double _Nagle_EN, double _Gordy_EN, double _Rahm_EN, double _Allen_EN, 
        double _Ghosh_EN, double _Mulliken_EN, double _AllredRochow_EN, 
        double _CotrellSutton_EN, double _Boeyens_VS_EN, double _electrophilicity,
        double _ionic_radius, double _crystal_radius, double _enthalpyAtomization, 
        double _specificHeat, double _critical_nuclear_charge, double _Tandon_EN, 
        double _LatticeConstant, double	_atomicRadius_Rahm)
    {
        symbol = _symbol;
        atomicNumber = _atomicNumber;
        atomicWeight = _atomicWeight;
        mendeleevNumber = _mendeleevNumber;
        atomicRadius = _atomicRadius;
        vdwRadius = _vdwRadius;
        
        LQuantumNumber = _LQuantumNumber;
        sOrbitalElectronNumber = _sOrbital;
        pOrbitalElectronNumber = _pOrbital;
        fOrbitalElectronNumber = _fOrbital;
        dOrbitalElectronNumber = _dOrbital;
        ionizationEnergy = _ionizationEnergy;
        dipolePolarizability = _dipolePolarizability;
        heatOfFusion = _heatOfFusion;
        heatOfVaporization = _heatOfVaporization;
        thermalConductivity = _thermalConductivity;
        cohesiveEnergy = _cohesiveEnergy;
        VEC = _VEC;
        MP = _MP;
        BP = _BP;
        density = _density;
        group = _group;
        electronAffinity = _electronAffinity;
        period = _period;
        spaceGroupNumber = _SGN;
        specificHeat = _specificHeat;
        latticeConstant = _LatticeConstant;
        atomicRadiusRahm = _atomicRadius_Rahm;
        critical_nuclear_charge = _critical_nuclear_charge;
        
        
        radiusZunger = _zungerRadius;
        bulkModulus = _bulkModulus;
        atomicVolume = _atomicVolume;
        covalent_radius = _covalent_radius;
        ionicRadius = _ionic_radius;
        crystalRadius = _crystal_radius;
        
        enthalpyAtomization = _enthalpyAtomization;
        electronDensity = _electronDensity;
        magneticMoment = _magneticMoment;
        workFunction = _workFunction;
        num_f_unfilled_valence_electrons = _num_f_unfilled_valence_electrons;
        num_p_unfilled_valence_electrons = _num_p_unfilled_valence_electrons;
        num_s_unfilled_valence_electrons = _num_s_unfilled_valence_electrons;
        num_d_unfilled_valence_electrons = _num_d_unfilled_valence_electrons;
        radius_d_orbital = _radius_d_orbital;
        radius_p_orbital = _radius_p_orbital;
        radius_s_orbital = _radius_s_orbital;
        total_unfilled_valence_electrons = _total_unfilled_valence_electrons;
        total_valence_electrons = _total_valence_electrons;
        
        en_Ghosh = _Ghosh_EN;
        en_Mulliken = _Mulliken_EN;
        en_Nagle = _Nagle_EN;
        en_Gordy = _Gordy_EN;
        en_Rahm = _Rahm_EN;
        en_Allen = _Allen_EN;
        en_AllredRochow = _AllredRochow_EN;
        en_CotrellSutton = _CotrellSutton_EN;
        en_Boeyens = _Boeyens_VS_EN;
        electrophilicity = _electrophilicity;
        en_MartynovBatsanov = _electronegativityMB;
        en_Pauling = _electronegativityPauling;
        en_Tandon = _Tandon_EN;
    }

    public double getAtomicVolume()
    {
        return atomicVolume;
    }

    public void setAtomicVolume(double atomicVolume)
    {
        this.atomicVolume = atomicVolume;
    }

    public double getBulkModulus()
    {
        return bulkModulus;
    }

    public void setBulkModulus(double bulkModulus)
    {
        this.bulkModulus = bulkModulus;
    }

    public int getSpaceGroupNumber()
    {
        return spaceGroupNumber;
    }

    public void setSpaceGroupNumber(int spaceGroupNumber)
    {
        this.spaceGroupNumber = spaceGroupNumber;
    }

    public double getEN_MartynovBatsanov()
    {
        return en_MartynovBatsanov;
    }

    public void setEN_MartynovBatsanov(double en_MartynovBatsanov)
    {
        this.en_MartynovBatsanov = en_MartynovBatsanov;
    }

    public double getRadiusZunger()
    {
        return radiusZunger;
    }

    public void setRadiusZunger(double radiusZunger)
    {
        this.radiusZunger = radiusZunger;
    }
    
    public String getSymbol()
    {
        return symbol;
    }

    public void setSymbol(String symbol)
    {
        this.symbol = symbol;
    }

    public int getAtomicNumber()
    {
        return atomicNumber;
    }

    public void setAtomicNumber(int atomicNumber)
    {
        this.atomicNumber = atomicNumber;
    }

    public double getAtomicWeight()
    {
        return atomicWeight;
    }

    public void setAtomicWeight(double atomicWeight)
    {
        this.atomicWeight = atomicWeight;
    }

    public double getVEC()
    {
        return VEC;
    }

    public void setVEC(double VEC)
    {
        this.VEC = VEC;
    }

    public double getAtomicRadius()
    {
        return atomicRadius;
    }

    public void setAtomicRadius(double atomicRadius)
    {
        this.atomicRadius = atomicRadius;
    }

    public double getEN_Pauling()
    {
        return en_Pauling;
    }

    public void setEN_Pauling(double en_Pauling)
    {
        this.en_Pauling = en_Pauling;
    }

    public double getMP()
    {
        return MP;
    }

    public void setMP(double MP)
    {
        this.MP = MP;
    }

    public double getBP()
    {
        return BP;
    }

    public void setBP(double BP)
    {
        this.BP = BP;
    }

    public double getThermalConductivity()
    {
        return thermalConductivity;
    }

    public void setThermalConductivity(double thermalConductivity)
    {
        this.thermalConductivity = thermalConductivity;
    }

    public double getCohesiveEnergy()
    {
        return cohesiveEnergy;
    }

    public void setCohesiveEnergy(double cohesiveEnergy)
    {
        this.cohesiveEnergy = cohesiveEnergy;
    }

    public double getIonizationEnergy()
    {
        return ionizationEnergy;
    }

    public void setIonizationEnergy(double ionizationEnergy)
    {
        this.ionizationEnergy = ionizationEnergy;
    }

    public double getDipolePolarizability()
    {
        return dipolePolarizability;
    }

    public void setDipolePolarizability(double dipolePolarizability)
    {
        this.dipolePolarizability = dipolePolarizability;
    }

    public double getHeatOfFusion()
    {
        return heatOfFusion;
    }

    public void setHeatOfFusion(double heatOfFusion)
    {
        this.heatOfFusion = heatOfFusion;
    }

    public double getHeatOfVaporization()
    {
        return heatOfVaporization;
    }

    public void setHeatOfVaporization(double heatOfVaporization)
    {
        this.heatOfVaporization = heatOfVaporization;
    }

    public int getMendeleevNumber()
    {
        return mendeleevNumber;
    }

    public void setMendeleevNumber(int mendeleevNumber)
    {
        this.mendeleevNumber = mendeleevNumber;
    }

    public int getLQuantumNumber()
    {
        return LQuantumNumber;
    }

    public void setLQuantumNumber(int LQuantumNumber)
    {
        this.LQuantumNumber = LQuantumNumber;
    }

    public double getDensity()
    {
        return density;
    }

    public void setDensity(double density)
    {
        this.density = density;
    }

    public double getElectronAffinity()
    {
        return electronAffinity;
    }

    public void setElectronAffinity(double electronAffinity)
    {
        this.electronAffinity = electronAffinity;
    }

    public int getsOrbitalElectronNumber()
    {
        return sOrbitalElectronNumber;
    }

    public void setsOrbitalElectronNumber(int sOrbitalElectronNumber)
    {
        this.sOrbitalElectronNumber = sOrbitalElectronNumber;
    }

    public int getdOrbitalElectronNumber()
    {
        return dOrbitalElectronNumber;
    }

    public void setdOrbitalElectronNumber(int dOrbitalElectronNumber)
    {
        this.dOrbitalElectronNumber = dOrbitalElectronNumber;
    }

    public int getfOrbitalElectronNumber()
    {
        return fOrbitalElectronNumber;
    }

    public void setfOrbitalElectronNumber(int fOrbitalElectronNumber)
    {
        this.fOrbitalElectronNumber = fOrbitalElectronNumber;
    }

    public int getpOrbitalElectronNumber()
    {
        return pOrbitalElectronNumber;
    }

    public void setpOrbitalElectronNumber(int pOrbitalElectronNumber)
    {
        this.pOrbitalElectronNumber = pOrbitalElectronNumber;
    }

    public double getVdwRadius()
    {
        return vdwRadius;
    }

    public void setVdwRadius(double vdwRadius)
    {
        this.vdwRadius = vdwRadius;
    }

    public int getValence()
    {
        return valence;
    }

    public void setValence(int valence)
    {
        this.valence = valence;
    }

    public double getGroup()
    {
        return group;
    }

    public void setGroup(double group)
    {
        this.group = group;
    }

    public double getPeriod()
    {
        return period;
    }

    public void setPeriod(double period)
    {
        this.period = period;
    }

    public double getCovalent_radius()
    {
        return covalent_radius;
    }

    public void setCovalent_radius(double covalent_radius)
    {
        this.covalent_radius = covalent_radius;
    }

    public double getElectronDensity()
    {
        return electronDensity;
    }

    public void setElectronDensity(double electronDensity)
    {
        this.electronDensity = electronDensity;
    }

    public double getMagneticMoment()
    {
        return magneticMoment;
    }

    public void setMagneticMoment(double magneticMoment)
    {
        this.magneticMoment = magneticMoment;
    }

    public double getWorkFunction()
    {
        return workFunction;
    }

    public void setWorkFunction(double workFunction)
    {
        this.workFunction = workFunction;
    }

    public int getNum_f_unfilled_valence_electrons()
    {
        return num_f_unfilled_valence_electrons;
    }

    public void setNum_f_unfilled_valence_electrons(
                                                    int num_f_unfilled_valence_electrons)
    {
        this.num_f_unfilled_valence_electrons = num_f_unfilled_valence_electrons;
    }

    public int getNum_p_unfilled_valence_electrons()
    {
        return num_p_unfilled_valence_electrons;
    }

    public void setNum_p_unfilled_valence_electrons(
                                                    int num_p_unfilled_valence_electrons)
    {
        this.num_p_unfilled_valence_electrons = num_p_unfilled_valence_electrons;
    }

    public int getNum_d_unfilled_valence_electrons()
    {
        return num_d_unfilled_valence_electrons;
    }

    public void setNum_d_unfilled_valence_electrons(
                                                    int num_d_unfilled_valence_electrons)
    {
        this.num_d_unfilled_valence_electrons = num_d_unfilled_valence_electrons;
    }
    
    

    public int getNum_s_unfilled_valence_electrons()
    {
        return num_s_unfilled_valence_electrons;
    }

    public void setNum_s_unfilled_valence_electrons(
                                                    int num_s_unfilled_valence_electrons)
    {
        this.num_s_unfilled_valence_electrons = num_s_unfilled_valence_electrons;
    }

    public double getRadius_d_orbital()
    {
        return radius_d_orbital;
    }

    public void setRadius_d_orbital(double radius_d_orbital)
    {
        this.radius_d_orbital = radius_d_orbital;
    }

    public double getRadius_p_orbital()
    {
        return radius_p_orbital;
    }

    public void setRadius_p_orbital(double radius_p_orbital)
    {
        this.radius_p_orbital = radius_p_orbital;
    }

    public double getRadius_s_orbital()
    {
        return radius_s_orbital;
    }

    public void setRadius_s_orbital(double radius_s_orbital)
    {
        this.radius_s_orbital = radius_s_orbital;
    }

    public int getTotal_unfilled_valence_electrons()
    {
        return total_unfilled_valence_electrons;
    }

    public void setTotal_unfilled_valence_electrons(
                                                    int total_unfilled_valence_electrons)
    {
        this.total_unfilled_valence_electrons = total_unfilled_valence_electrons;
    }

    public int getTotal_valence_electrons()
    {
        return total_valence_electrons;
    }

    public void setTotal_valence_electrons(int total_valence_electrons)
    {
        this.total_valence_electrons = total_valence_electrons;
    }

    public double getEN_Ghosh()
    {
        return en_Ghosh;
    }

    public void setEN_Ghosh(double en_Ghosh)
    {
        this.en_Ghosh = en_Ghosh;
    }

    public double getEN_Mulliken()
    {
        return en_Mulliken;
    }

    public void setEN_Mulliken(double en_Mulliken)
    {
        this.en_Mulliken = en_Mulliken;
    }

    public double getEN_Nagle()
    {
        return en_Nagle;
    }

    public void setEN_Nagle(double en_Nagle)
    {
        this.en_Nagle = en_Nagle;
    }

    public double getEN_Gordy()
    {
        return en_Gordy;
    }

    public void setEN_Gordy(double en_Gordy)
    {
        this.en_Gordy = en_Gordy;
    }

    public double getEN_Rahm()
    {
        return en_Rahm;
    }

    public void setEN_Rahm(double en_Rahm)
    {
        this.en_Rahm = en_Rahm;
    }

    public double getEN_Allen()
    {
        return en_Allen;
    }

    public void setEN_Allen(double en_Allen)
    {
        this.en_Allen = en_Allen;
    }

    public double getEN_AllredRochow()
    {
        return en_AllredRochow;
    }

    public void setEN_AllredRochow(double en_AllredRochow)
    {
        this.en_AllredRochow = en_AllredRochow;
    }

    public double getEN_CotrellSutton()
    {
        return en_CotrellSutton;
    }

    public void setEN_CotrellSutton(double en_CotrellSutton)
    {
        this.en_CotrellSutton = en_CotrellSutton;
    }

    public double getEN_Boeyens()
    {
        return en_Boeyens;
    }

    public void setEN_Boeyens(double en_Boeyens)
    {
        this.en_Boeyens = en_Boeyens;
    }

    public double getElectrophilicity()
    {
        return electrophilicity;
    }

    public void setElectrophilicity(double electrophilicity)
    {
        this.electrophilicity = electrophilicity;
    }

    public double getIonicRadius()
    {
        return ionicRadius;
    }

    public void setIonicRadius(double ionicRadius)
    {
        this.ionicRadius = ionicRadius;
    }

    public double getCrystalRadius()
    {
        return crystalRadius;
    }

    public void setCrystalRadius(double crystalRadius)
    {
        this.crystalRadius = crystalRadius;
    }

    public double getEnthalpyAtomization()
    {
        return enthalpyAtomization;
    }

    public void setEnthalpyAtomization(double enthalpyAtomization)
    {
        this.enthalpyAtomization = enthalpyAtomization;
    }

    public double getSpecificHeat()
    {
        return specificHeat;
    }

    public void setSpecificHeat(double specificHeat)
    {
        this.specificHeat = specificHeat;
    }

    public double getCriticalNuclearCharge()
    {
        return critical_nuclear_charge;
    }

    public void setCriticalNuclearCharge(double critical_nuclear_charge)
    {
        this.critical_nuclear_charge = critical_nuclear_charge;
    }

    public double getEN_Tandon()
    {
        return en_Tandon;
    }

    public void setEN_Tandon(double en_Tandon)
    {
        this.en_Tandon = en_Tandon;
    }

    public double getLatticeConstant()
    {
        return latticeConstant;
    }

    public void setLatticeConstant(double latticeConstant)
    {
        this.latticeConstant = latticeConstant;
    }

    public double getAtomicRadiusRahm()
    {
        return atomicRadiusRahm;
    }

    public void setAtomicRadiusRahm(double atomicRadiusRahm)
    {
        this.atomicRadiusRahm = atomicRadiusRahm;
    }
    
    
    
}

