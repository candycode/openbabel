/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 1.3.31
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */


public class OBMolAtomBFSIter {
  private long swigCPtr;
  protected boolean swigCMemOwn;

  protected OBMolAtomBFSIter(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(OBMolAtomBFSIter obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if(swigCPtr != 0 && swigCMemOwn) {
      swigCMemOwn = false;
      openbabelJNI.delete_OBMolAtomBFSIter(swigCPtr);
    }
    swigCPtr = 0;
  }

  public OBMolAtomBFSIter() {
    this(openbabelJNI.new_OBMolAtomBFSIter__SWIG_0(), true);
  }

  public OBMolAtomBFSIter(OBMol mol) {
    this(openbabelJNI.new_OBMolAtomBFSIter__SWIG_1(OBMol.getCPtr(mol), mol), true);
  }

  public OBMolAtomBFSIter(OBMolAtomBFSIter ai) {
    this(openbabelJNI.new_OBMolAtomBFSIter__SWIG_2(OBMolAtomBFSIter.getCPtr(ai), ai), true);
  }

  public OBAtom __deref__() {
    long cPtr = openbabelJNI.OBMolAtomBFSIter___deref__(swigCPtr, this);
    return (cPtr == 0) ? null : new OBAtom(cPtr, false);
  }

  public OBAtom __ref__() {
    return new OBAtom(openbabelJNI.OBMolAtomBFSIter___ref__(swigCPtr, this), false);
  }

  public void setVisit(boolean value) {
    openbabelJNI.OBMolAtomBFSIter_Visit_set(swigCPtr, this, value);
  }

  public boolean getVisit() {
    return openbabelJNI.OBMolAtomBFSIter_Visit_get(swigCPtr, this);
  }

  public boolean Clear() {
    return openbabelJNI.OBMolAtomBFSIter_Clear(swigCPtr, this);
  }

  public void SetIdx(int idx) {
    openbabelJNI.OBMolAtomBFSIter_SetIdx(swigCPtr, this, idx);
  }

  public void SetHyb(int hyb) {
    openbabelJNI.OBMolAtomBFSIter_SetHyb(swigCPtr, this, hyb);
  }

  public void SetAtomicNum(int atomicnum) {
    openbabelJNI.OBMolAtomBFSIter_SetAtomicNum(swigCPtr, this, atomicnum);
  }

  public void SetIsotope(long iso) {
    openbabelJNI.OBMolAtomBFSIter_SetIsotope(swigCPtr, this, iso);
  }

  public void SetImplicitValence(int val) {
    openbabelJNI.OBMolAtomBFSIter_SetImplicitValence(swigCPtr, this, val);
  }

  public void IncrementImplicitValence() {
    openbabelJNI.OBMolAtomBFSIter_IncrementImplicitValence(swigCPtr, this);
  }

  public void DecrementImplicitValence() {
    openbabelJNI.OBMolAtomBFSIter_DecrementImplicitValence(swigCPtr, this);
  }

  public void SetFormalCharge(int fcharge) {
    openbabelJNI.OBMolAtomBFSIter_SetFormalCharge(swigCPtr, this, fcharge);
  }

  public void SetSpinMultiplicity(short spin) {
    openbabelJNI.OBMolAtomBFSIter_SetSpinMultiplicity(swigCPtr, this, spin);
  }

  public void SetType(String type) {
    openbabelJNI.OBMolAtomBFSIter_SetType__SWIG_0(swigCPtr, this, type);
  }

  public void SetType(SWIGTYPE_p_std__string type) {
    openbabelJNI.OBMolAtomBFSIter_SetType__SWIG_1(swigCPtr, this, SWIGTYPE_p_std__string.getCPtr(type));
  }

  public void SetPartialCharge(double pcharge) {
    openbabelJNI.OBMolAtomBFSIter_SetPartialCharge(swigCPtr, this, pcharge);
  }

  public void SetVector(vector3 v) {
    openbabelJNI.OBMolAtomBFSIter_SetVector__SWIG_0(swigCPtr, this, vector3.getCPtr(v), v);
  }

  public void SetVector(double x, double y, double z) {
    openbabelJNI.OBMolAtomBFSIter_SetVector__SWIG_1(swigCPtr, this, x, y, z);
  }

  public void SetVector() {
    openbabelJNI.OBMolAtomBFSIter_SetVector__SWIG_2(swigCPtr, this);
  }

  public void SetCoordPtr(SWIGTYPE_p_p_double c) {
    openbabelJNI.OBMolAtomBFSIter_SetCoordPtr(swigCPtr, this, SWIGTYPE_p_p_double.getCPtr(c));
  }

  public void SetResidue(OBResidue res) {
    openbabelJNI.OBMolAtomBFSIter_SetResidue(swigCPtr, this, OBResidue.getCPtr(res), res);
  }

  public void SetParent(OBMol ptr) {
    openbabelJNI.OBMolAtomBFSIter_SetParent(swigCPtr, this, OBMol.getCPtr(ptr), ptr);
  }

  public void SetAromatic() {
    openbabelJNI.OBMolAtomBFSIter_SetAromatic(swigCPtr, this);
  }

  public void UnsetAromatic() {
    openbabelJNI.OBMolAtomBFSIter_UnsetAromatic(swigCPtr, this);
  }

  public void SetClockwiseStereo() {
    openbabelJNI.OBMolAtomBFSIter_SetClockwiseStereo(swigCPtr, this);
  }

  public void SetAntiClockwiseStereo() {
    openbabelJNI.OBMolAtomBFSIter_SetAntiClockwiseStereo(swigCPtr, this);
  }

  public void SetPositiveStereo() {
    openbabelJNI.OBMolAtomBFSIter_SetPositiveStereo(swigCPtr, this);
  }

  public void SetNegativeStereo() {
    openbabelJNI.OBMolAtomBFSIter_SetNegativeStereo(swigCPtr, this);
  }

  public void UnsetStereo() {
    openbabelJNI.OBMolAtomBFSIter_UnsetStereo(swigCPtr, this);
  }

  public void SetInRing() {
    openbabelJNI.OBMolAtomBFSIter_SetInRing(swigCPtr, this);
  }

  public void SetChiral() {
    openbabelJNI.OBMolAtomBFSIter_SetChiral(swigCPtr, this);
  }

  public void ClearCoordPtr() {
    openbabelJNI.OBMolAtomBFSIter_ClearCoordPtr(swigCPtr, this);
  }

  public int GetFormalCharge() {
    return openbabelJNI.OBMolAtomBFSIter_GetFormalCharge(swigCPtr, this);
  }

  public long GetAtomicNum() {
    return openbabelJNI.OBMolAtomBFSIter_GetAtomicNum(swigCPtr, this);
  }

  public int GetIsotope() {
    return openbabelJNI.OBMolAtomBFSIter_GetIsotope(swigCPtr, this);
  }

  public int GetSpinMultiplicity() {
    return openbabelJNI.OBMolAtomBFSIter_GetSpinMultiplicity(swigCPtr, this);
  }

  public double GetAtomicMass() {
    return openbabelJNI.OBMolAtomBFSIter_GetAtomicMass(swigCPtr, this);
  }

  public double GetExactMass() {
    return openbabelJNI.OBMolAtomBFSIter_GetExactMass(swigCPtr, this);
  }

  public long GetIdx() {
    return openbabelJNI.OBMolAtomBFSIter_GetIdx(swigCPtr, this);
  }

  public long GetCoordinateIdx() {
    return openbabelJNI.OBMolAtomBFSIter_GetCoordinateIdx(swigCPtr, this);
  }

  public long GetCIdx() {
    return openbabelJNI.OBMolAtomBFSIter_GetCIdx(swigCPtr, this);
  }

  public long GetValence() {
    return openbabelJNI.OBMolAtomBFSIter_GetValence(swigCPtr, this);
  }

  public long GetHyb() {
    return openbabelJNI.OBMolAtomBFSIter_GetHyb(swigCPtr, this);
  }

  public long GetImplicitValence() {
    return openbabelJNI.OBMolAtomBFSIter_GetImplicitValence(swigCPtr, this);
  }

  public long GetHvyValence() {
    return openbabelJNI.OBMolAtomBFSIter_GetHvyValence(swigCPtr, this);
  }

  public long GetHeteroValence() {
    return openbabelJNI.OBMolAtomBFSIter_GetHeteroValence(swigCPtr, this);
  }

  public String GetType() {
    return openbabelJNI.OBMolAtomBFSIter_GetType(swigCPtr, this);
  }

  public double GetX() {
    return openbabelJNI.OBMolAtomBFSIter_GetX(swigCPtr, this);
  }

  public double GetY() {
    return openbabelJNI.OBMolAtomBFSIter_GetY(swigCPtr, this);
  }

  public double GetZ() {
    return openbabelJNI.OBMolAtomBFSIter_GetZ(swigCPtr, this);
  }

  public double x() {
    return openbabelJNI.OBMolAtomBFSIter_x(swigCPtr, this);
  }

  public double y() {
    return openbabelJNI.OBMolAtomBFSIter_y(swigCPtr, this);
  }

  public double z() {
    return openbabelJNI.OBMolAtomBFSIter_z(swigCPtr, this);
  }

  public SWIGTYPE_p_double GetCoordinate() {
    long cPtr = openbabelJNI.OBMolAtomBFSIter_GetCoordinate(swigCPtr, this);
    return (cPtr == 0) ? null : new SWIGTYPE_p_double(cPtr, false);
  }

  public vector3 GetVector() {
    return new vector3(openbabelJNI.OBMolAtomBFSIter_GetVector__SWIG_0(swigCPtr, this), false);
  }

  public double GetPartialCharge() {
    return openbabelJNI.OBMolAtomBFSIter_GetPartialCharge(swigCPtr, this);
  }

  public OBResidue GetResidue() {
    long cPtr = openbabelJNI.OBMolAtomBFSIter_GetResidue(swigCPtr, this);
    return (cPtr == 0) ? null : new OBResidue(cPtr, false);
  }

  public OBMol GetParent() {
    long cPtr = openbabelJNI.OBMolAtomBFSIter_GetParent(swigCPtr, this);
    return (cPtr == 0) ? null : new OBMol(cPtr, false);
  }

  public boolean GetNewBondVector(vector3 v, double length) {
    return openbabelJNI.OBMolAtomBFSIter_GetNewBondVector(swigCPtr, this, vector3.getCPtr(v), v, length);
  }

  public OBBond GetBond(OBAtom arg0) {
    long cPtr = openbabelJNI.OBMolAtomBFSIter_GetBond(swigCPtr, this, OBAtom.getCPtr(arg0), arg0);
    return (cPtr == 0) ? null : new OBBond(cPtr, false);
  }

  public OBAtom GetNextAtom() {
    long cPtr = openbabelJNI.OBMolAtomBFSIter_GetNextAtom(swigCPtr, this);
    return (cPtr == 0) ? null : new OBAtom(cPtr, false);
  }

  public SWIGTYPE_p_std__vectorTOpenBabel__OBBond_p_t__iterator BeginBonds() {
    return new SWIGTYPE_p_std__vectorTOpenBabel__OBBond_p_t__iterator(openbabelJNI.OBMolAtomBFSIter_BeginBonds(swigCPtr, this), true);
  }

  public SWIGTYPE_p_std__vectorTOpenBabel__OBBond_p_t__iterator EndBonds() {
    return new SWIGTYPE_p_std__vectorTOpenBabel__OBBond_p_t__iterator(openbabelJNI.OBMolAtomBFSIter_EndBonds(swigCPtr, this), true);
  }

  public OBBond BeginBond(SWIGTYPE_p_std__vectorTOpenBabel__OBBond_p_t__iterator i) {
    long cPtr = openbabelJNI.OBMolAtomBFSIter_BeginBond(swigCPtr, this, SWIGTYPE_p_std__vectorTOpenBabel__OBBond_p_t__iterator.getCPtr(i));
    return (cPtr == 0) ? null : new OBBond(cPtr, false);
  }

  public OBBond NextBond(SWIGTYPE_p_std__vectorTOpenBabel__OBBond_p_t__iterator i) {
    long cPtr = openbabelJNI.OBMolAtomBFSIter_NextBond(swigCPtr, this, SWIGTYPE_p_std__vectorTOpenBabel__OBBond_p_t__iterator.getCPtr(i));
    return (cPtr == 0) ? null : new OBBond(cPtr, false);
  }

  public OBAtom BeginNbrAtom(SWIGTYPE_p_std__vectorTOpenBabel__OBBond_p_t__iterator i) {
    long cPtr = openbabelJNI.OBMolAtomBFSIter_BeginNbrAtom(swigCPtr, this, SWIGTYPE_p_std__vectorTOpenBabel__OBBond_p_t__iterator.getCPtr(i));
    return (cPtr == 0) ? null : new OBAtom(cPtr, false);
  }

  public OBAtom NextNbrAtom(SWIGTYPE_p_std__vectorTOpenBabel__OBBond_p_t__iterator i) {
    long cPtr = openbabelJNI.OBMolAtomBFSIter_NextNbrAtom(swigCPtr, this, SWIGTYPE_p_std__vectorTOpenBabel__OBBond_p_t__iterator.getCPtr(i));
    return (cPtr == 0) ? null : new OBAtom(cPtr, false);
  }

  public double GetDistance(int index) {
    return openbabelJNI.OBMolAtomBFSIter_GetDistance__SWIG_0(swigCPtr, this, index);
  }

  public double GetDistance(OBAtom arg0) {
    return openbabelJNI.OBMolAtomBFSIter_GetDistance__SWIG_1(swigCPtr, this, OBAtom.getCPtr(arg0), arg0);
  }

  public double GetAngle(int b, int c) {
    return openbabelJNI.OBMolAtomBFSIter_GetAngle__SWIG_0(swigCPtr, this, b, c);
  }

  public double GetAngle(OBAtom b, OBAtom c) {
    return openbabelJNI.OBMolAtomBFSIter_GetAngle__SWIG_1(swigCPtr, this, OBAtom.getCPtr(b), b, OBAtom.getCPtr(c), c);
  }

  public void NewResidue() {
    openbabelJNI.OBMolAtomBFSIter_NewResidue(swigCPtr, this);
  }

  public void AddResidue(OBResidue res) {
    openbabelJNI.OBMolAtomBFSIter_AddResidue(swigCPtr, this, OBResidue.getCPtr(res), res);
  }

  public void DeleteResidue() {
    openbabelJNI.OBMolAtomBFSIter_DeleteResidue(swigCPtr, this);
  }

  public void AddBond(OBBond bond) {
    openbabelJNI.OBMolAtomBFSIter_AddBond(swigCPtr, this, OBBond.getCPtr(bond), bond);
  }

  public void InsertBond(SWIGTYPE_p_std__vectorTOpenBabel__OBBond_p_t__iterator i, OBBond bond) {
    openbabelJNI.OBMolAtomBFSIter_InsertBond(swigCPtr, this, SWIGTYPE_p_std__vectorTOpenBabel__OBBond_p_t__iterator.getCPtr(i), OBBond.getCPtr(bond), bond);
  }

  public boolean DeleteBond(OBBond bond) {
    return openbabelJNI.OBMolAtomBFSIter_DeleteBond(swigCPtr, this, OBBond.getCPtr(bond), bond);
  }

  public void ClearBond() {
    openbabelJNI.OBMolAtomBFSIter_ClearBond(swigCPtr, this);
  }

  public boolean HtoMethyl() {
    return openbabelJNI.OBMolAtomBFSIter_HtoMethyl(swigCPtr, this);
  }

  public boolean SetHybAndGeom(int arg0) {
    return openbabelJNI.OBMolAtomBFSIter_SetHybAndGeom(swigCPtr, this, arg0);
  }

  public void ForceNoH() {
    openbabelJNI.OBMolAtomBFSIter_ForceNoH(swigCPtr, this);
  }

  public boolean HasNoHForced() {
    return openbabelJNI.OBMolAtomBFSIter_HasNoHForced(swigCPtr, this);
  }

  public long CountFreeOxygens() {
    return openbabelJNI.OBMolAtomBFSIter_CountFreeOxygens(swigCPtr, this);
  }

  public long ImplicitHydrogenCount() {
    return openbabelJNI.OBMolAtomBFSIter_ImplicitHydrogenCount(swigCPtr, this);
  }

  public long ExplicitHydrogenCount(boolean ExcludeIsotopes) {
    return openbabelJNI.OBMolAtomBFSIter_ExplicitHydrogenCount__SWIG_0(swigCPtr, this, ExcludeIsotopes);
  }

  public long ExplicitHydrogenCount() {
    return openbabelJNI.OBMolAtomBFSIter_ExplicitHydrogenCount__SWIG_1(swigCPtr, this);
  }

  public long MemberOfRingCount() {
    return openbabelJNI.OBMolAtomBFSIter_MemberOfRingCount(swigCPtr, this);
  }

  public long MemberOfRingSize() {
    return openbabelJNI.OBMolAtomBFSIter_MemberOfRingSize(swigCPtr, this);
  }

  public long CountRingBonds() {
    return openbabelJNI.OBMolAtomBFSIter_CountRingBonds(swigCPtr, this);
  }

  public double SmallestBondAngle() {
    return openbabelJNI.OBMolAtomBFSIter_SmallestBondAngle(swigCPtr, this);
  }

  public double AverageBondAngle() {
    return openbabelJNI.OBMolAtomBFSIter_AverageBondAngle(swigCPtr, this);
  }

  public long BOSum() {
    return openbabelJNI.OBMolAtomBFSIter_BOSum(swigCPtr, this);
  }

  public long KBOSum() {
    return openbabelJNI.OBMolAtomBFSIter_KBOSum(swigCPtr, this);
  }

  public boolean HasResidue() {
    return openbabelJNI.OBMolAtomBFSIter_HasResidue(swigCPtr, this);
  }

  public boolean IsHydrogen() {
    return openbabelJNI.OBMolAtomBFSIter_IsHydrogen(swigCPtr, this);
  }

  public boolean IsCarbon() {
    return openbabelJNI.OBMolAtomBFSIter_IsCarbon(swigCPtr, this);
  }

  public boolean IsNitrogen() {
    return openbabelJNI.OBMolAtomBFSIter_IsNitrogen(swigCPtr, this);
  }

  public boolean IsOxygen() {
    return openbabelJNI.OBMolAtomBFSIter_IsOxygen(swigCPtr, this);
  }

  public boolean IsSulfur() {
    return openbabelJNI.OBMolAtomBFSIter_IsSulfur(swigCPtr, this);
  }

  public boolean IsPhosphorus() {
    return openbabelJNI.OBMolAtomBFSIter_IsPhosphorus(swigCPtr, this);
  }

  public boolean IsAromatic() {
    return openbabelJNI.OBMolAtomBFSIter_IsAromatic(swigCPtr, this);
  }

  public boolean IsInRing() {
    return openbabelJNI.OBMolAtomBFSIter_IsInRing(swigCPtr, this);
  }

  public boolean IsInRingSize(int arg0) {
    return openbabelJNI.OBMolAtomBFSIter_IsInRingSize(swigCPtr, this, arg0);
  }

  public boolean IsHeteroatom() {
    return openbabelJNI.OBMolAtomBFSIter_IsHeteroatom(swigCPtr, this);
  }

  public boolean IsNotCorH() {
    return openbabelJNI.OBMolAtomBFSIter_IsNotCorH(swigCPtr, this);
  }

  public boolean IsConnected(OBAtom arg0) {
    return openbabelJNI.OBMolAtomBFSIter_IsConnected(swigCPtr, this, OBAtom.getCPtr(arg0), arg0);
  }

  public boolean IsOneThree(OBAtom arg0) {
    return openbabelJNI.OBMolAtomBFSIter_IsOneThree(swigCPtr, this, OBAtom.getCPtr(arg0), arg0);
  }

  public boolean IsOneFour(OBAtom arg0) {
    return openbabelJNI.OBMolAtomBFSIter_IsOneFour(swigCPtr, this, OBAtom.getCPtr(arg0), arg0);
  }

  public boolean IsCarboxylOxygen() {
    return openbabelJNI.OBMolAtomBFSIter_IsCarboxylOxygen(swigCPtr, this);
  }

  public boolean IsPhosphateOxygen() {
    return openbabelJNI.OBMolAtomBFSIter_IsPhosphateOxygen(swigCPtr, this);
  }

  public boolean IsSulfateOxygen() {
    return openbabelJNI.OBMolAtomBFSIter_IsSulfateOxygen(swigCPtr, this);
  }

  public boolean IsNitroOxygen() {
    return openbabelJNI.OBMolAtomBFSIter_IsNitroOxygen(swigCPtr, this);
  }

  public boolean IsAmideNitrogen() {
    return openbabelJNI.OBMolAtomBFSIter_IsAmideNitrogen(swigCPtr, this);
  }

  public boolean IsPolarHydrogen() {
    return openbabelJNI.OBMolAtomBFSIter_IsPolarHydrogen(swigCPtr, this);
  }

  public boolean IsNonPolarHydrogen() {
    return openbabelJNI.OBMolAtomBFSIter_IsNonPolarHydrogen(swigCPtr, this);
  }

  public boolean IsAromaticNOxide() {
    return openbabelJNI.OBMolAtomBFSIter_IsAromaticNOxide(swigCPtr, this);
  }

  public boolean IsChiral() {
    return openbabelJNI.OBMolAtomBFSIter_IsChiral(swigCPtr, this);
  }

  public boolean IsAxial() {
    return openbabelJNI.OBMolAtomBFSIter_IsAxial(swigCPtr, this);
  }

  public boolean IsClockwise() {
    return openbabelJNI.OBMolAtomBFSIter_IsClockwise(swigCPtr, this);
  }

  public boolean IsAntiClockwise() {
    return openbabelJNI.OBMolAtomBFSIter_IsAntiClockwise(swigCPtr, this);
  }

  public boolean IsPositiveStereo() {
    return openbabelJNI.OBMolAtomBFSIter_IsPositiveStereo(swigCPtr, this);
  }

  public boolean IsNegativeStereo() {
    return openbabelJNI.OBMolAtomBFSIter_IsNegativeStereo(swigCPtr, this);
  }

  public boolean HasChiralitySpecified() {
    return openbabelJNI.OBMolAtomBFSIter_HasChiralitySpecified(swigCPtr, this);
  }

  public boolean HasChiralVolume() {
    return openbabelJNI.OBMolAtomBFSIter_HasChiralVolume(swigCPtr, this);
  }

  public boolean IsHbondAcceptor() {
    return openbabelJNI.OBMolAtomBFSIter_IsHbondAcceptor(swigCPtr, this);
  }

  public boolean IsHbondDonor() {
    return openbabelJNI.OBMolAtomBFSIter_IsHbondDonor(swigCPtr, this);
  }

  public boolean IsHbondDonorH() {
    return openbabelJNI.OBMolAtomBFSIter_IsHbondDonorH(swigCPtr, this);
  }

  public boolean HasAlphaBetaUnsat(boolean includePandS) {
    return openbabelJNI.OBMolAtomBFSIter_HasAlphaBetaUnsat__SWIG_0(swigCPtr, this, includePandS);
  }

  public boolean HasAlphaBetaUnsat() {
    return openbabelJNI.OBMolAtomBFSIter_HasAlphaBetaUnsat__SWIG_1(swigCPtr, this);
  }

  public boolean HasBondOfOrder(long bo) {
    return openbabelJNI.OBMolAtomBFSIter_HasBondOfOrder(swigCPtr, this, bo);
  }

  public int CountBondsOfOrder(long bo) {
    return openbabelJNI.OBMolAtomBFSIter_CountBondsOfOrder(swigCPtr, this, bo);
  }

  public boolean HasNonSingleBond() {
    return openbabelJNI.OBMolAtomBFSIter_HasNonSingleBond(swigCPtr, this);
  }

  public boolean HasSingleBond() {
    return openbabelJNI.OBMolAtomBFSIter_HasSingleBond(swigCPtr, this);
  }

  public boolean HasDoubleBond() {
    return openbabelJNI.OBMolAtomBFSIter_HasDoubleBond(swigCPtr, this);
  }

  public boolean HasAromaticBond() {
    return openbabelJNI.OBMolAtomBFSIter_HasAromaticBond(swigCPtr, this);
  }

  public boolean MatchesSMARTS(String arg0) {
    return openbabelJNI.OBMolAtomBFSIter_MatchesSMARTS(swigCPtr, this, arg0);
  }

  public OBBase DoTransformations(SWIGTYPE_p_std__mapTstd__string_std__string_t arg0) {
    long cPtr = openbabelJNI.OBMolAtomBFSIter_DoTransformations(swigCPtr, this, SWIGTYPE_p_std__mapTstd__string_std__string_t.getCPtr(arg0));
    return (cPtr == 0) ? null : new OBBase(cPtr, false);
  }

  public String ClassDescription() {
    return openbabelJNI.OBMolAtomBFSIter_ClassDescription(swigCPtr, this);
  }

  public boolean HasData(long type) {
    return openbabelJNI.OBMolAtomBFSIter_HasData__SWIG_2(swigCPtr, this, type);
  }

  public void DeleteData(long type) {
    openbabelJNI.OBMolAtomBFSIter_DeleteData__SWIG_0(swigCPtr, this, type);
  }

  public void DeleteData(OBGenericData arg0) {
    openbabelJNI.OBMolAtomBFSIter_DeleteData__SWIG_1(swigCPtr, this, OBGenericData.getCPtr(arg0), arg0);
  }

  public void DeleteData(vectorData arg0) {
    openbabelJNI.OBMolAtomBFSIter_DeleteData__SWIG_2(swigCPtr, this, vectorData.getCPtr(arg0), arg0);
  }

  public boolean DeleteData(String s) {
    return openbabelJNI.OBMolAtomBFSIter_DeleteData__SWIG_3(swigCPtr, this, s);
  }

  public void SetData(OBGenericData d) {
    openbabelJNI.OBMolAtomBFSIter_SetData(swigCPtr, this, OBGenericData.getCPtr(d), d);
  }

  public long DataSize() {
    return openbabelJNI.OBMolAtomBFSIter_DataSize(swigCPtr, this);
  }

  public OBGenericData GetData(long type) {
    long cPtr = openbabelJNI.OBMolAtomBFSIter_GetData__SWIG_0(swigCPtr, this, type);
    return (cPtr == 0) ? null : new OBGenericData(cPtr, false);
  }

  public vectorData GetData() {
    return new vectorData(openbabelJNI.OBMolAtomBFSIter_GetData__SWIG_3(swigCPtr, this), false);
  }

  public vectorData GetData(DataOrigin source) {
    return new vectorData(openbabelJNI.OBMolAtomBFSIter_GetData__SWIG_4(swigCPtr, this, source.swigValue()), true);
  }

  public SWIGTYPE_p_std__vectorTOpenBabel__OBGenericData_p_t__iterator BeginData() {
    return new SWIGTYPE_p_std__vectorTOpenBabel__OBGenericData_p_t__iterator(openbabelJNI.OBMolAtomBFSIter_BeginData(swigCPtr, this), true);
  }

  public SWIGTYPE_p_std__vectorTOpenBabel__OBGenericData_p_t__iterator EndData() {
    return new SWIGTYPE_p_std__vectorTOpenBabel__OBGenericData_p_t__iterator(openbabelJNI.OBMolAtomBFSIter_EndData(swigCPtr, this), true);
  }

}
