/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 1.3.31
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */


public class OBResidue extends OBBase {
  private long swigCPtr;

  protected OBResidue(long cPtr, boolean cMemoryOwn) {
    super(openbabelJNI.SWIGOBResidueUpcast(cPtr), cMemoryOwn);
    swigCPtr = cPtr;
  }

  protected static long getCPtr(OBResidue obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if(swigCPtr != 0 && swigCMemOwn) {
      swigCMemOwn = false;
      openbabelJNI.delete_OBResidue(swigCPtr);
    }
    swigCPtr = 0;
    super.delete();
  }

  public OBResidue() {
    this(openbabelJNI.new_OBResidue__SWIG_0(), true);
  }

  public OBResidue(OBResidue arg0) {
    this(openbabelJNI.new_OBResidue__SWIG_1(OBResidue.getCPtr(arg0), arg0), true);
  }

  public void AddAtom(OBAtom atom) {
    openbabelJNI.OBResidue_AddAtom(swigCPtr, this, OBAtom.getCPtr(atom), atom);
  }

  public void InsertAtom(OBAtom atom) {
    openbabelJNI.OBResidue_InsertAtom(swigCPtr, this, OBAtom.getCPtr(atom), atom);
  }

  public void RemoveAtom(OBAtom atom) {
    openbabelJNI.OBResidue_RemoveAtom(swigCPtr, this, OBAtom.getCPtr(atom), atom);
  }

  public boolean Clear() {
    return openbabelJNI.OBResidue_Clear(swigCPtr, this);
  }

  public void SetName(String resname) {
    openbabelJNI.OBResidue_SetName(swigCPtr, this, resname);
  }

  public void SetNum(long resnum) {
    openbabelJNI.OBResidue_SetNum(swigCPtr, this, resnum);
  }

  public void SetChain(char chain) {
    openbabelJNI.OBResidue_SetChain(swigCPtr, this, chain);
  }

  public void SetChainNum(long chainnum) {
    openbabelJNI.OBResidue_SetChainNum(swigCPtr, this, chainnum);
  }

  public void SetIdx(long idx) {
    openbabelJNI.OBResidue_SetIdx(swigCPtr, this, idx);
  }

  public void SetAtomID(OBAtom atom, String id) {
    openbabelJNI.OBResidue_SetAtomID(swigCPtr, this, OBAtom.getCPtr(atom), atom, id);
  }

  public void SetHetAtom(OBAtom atom, boolean hetatm) {
    openbabelJNI.OBResidue_SetHetAtom(swigCPtr, this, OBAtom.getCPtr(atom), atom, hetatm);
  }

  public void SetSerialNum(OBAtom atom, long sernum) {
    openbabelJNI.OBResidue_SetSerialNum(swigCPtr, this, OBAtom.getCPtr(atom), atom, sernum);
  }

  public String GetName() {
    return openbabelJNI.OBResidue_GetName(swigCPtr, this);
  }

  public long GetNum() {
    return openbabelJNI.OBResidue_GetNum(swigCPtr, this);
  }

  public long GetNumAtoms() {
    return openbabelJNI.OBResidue_GetNumAtoms(swigCPtr, this);
  }

  public char GetChain() {
    return openbabelJNI.OBResidue_GetChain(swigCPtr, this);
  }

  public long GetChainNum() {
    return openbabelJNI.OBResidue_GetChainNum(swigCPtr, this);
  }

  public long GetIdx() {
    return openbabelJNI.OBResidue_GetIdx(swigCPtr, this);
  }

  public long GetResKey() {
    return openbabelJNI.OBResidue_GetResKey(swigCPtr, this);
  }

  public SWIGTYPE_p_std__vectorTOpenBabel__OBAtom_p_t GetAtoms() {
    return new SWIGTYPE_p_std__vectorTOpenBabel__OBAtom_p_t(openbabelJNI.OBResidue_GetAtoms(swigCPtr, this), true);
  }

  public SWIGTYPE_p_std__vectorTOpenBabel__OBBond_p_t GetBonds(boolean exterior) {
    return new SWIGTYPE_p_std__vectorTOpenBabel__OBBond_p_t(openbabelJNI.OBResidue_GetBonds__SWIG_0(swigCPtr, this, exterior), true);
  }

  public SWIGTYPE_p_std__vectorTOpenBabel__OBBond_p_t GetBonds() {
    return new SWIGTYPE_p_std__vectorTOpenBabel__OBBond_p_t(openbabelJNI.OBResidue_GetBonds__SWIG_1(swigCPtr, this), true);
  }

  public String GetAtomID(OBAtom atom) {
    return openbabelJNI.OBResidue_GetAtomID(swigCPtr, this, OBAtom.getCPtr(atom), atom);
  }

  public long GetSerialNum(OBAtom atom) {
    return openbabelJNI.OBResidue_GetSerialNum(swigCPtr, this, OBAtom.getCPtr(atom), atom);
  }

  public boolean GetAminoAcidProperty(int arg0) {
    return openbabelJNI.OBResidue_GetAminoAcidProperty(swigCPtr, this, arg0);
  }

  public boolean GetAtomProperty(OBAtom a, int arg1) {
    return openbabelJNI.OBResidue_GetAtomProperty(swigCPtr, this, OBAtom.getCPtr(a), a, arg1);
  }

  public boolean GetResidueProperty(int arg0) {
    return openbabelJNI.OBResidue_GetResidueProperty(swigCPtr, this, arg0);
  }

  public boolean IsHetAtom(OBAtom atom) {
    return openbabelJNI.OBResidue_IsHetAtom(swigCPtr, this, OBAtom.getCPtr(atom), atom);
  }

  public boolean IsResidueType(int arg0) {
    return openbabelJNI.OBResidue_IsResidueType(swigCPtr, this, arg0);
  }

  public SWIGTYPE_p_std__vectorTOpenBabel__OBAtom_p_t__iterator BeginAtoms() {
    return new SWIGTYPE_p_std__vectorTOpenBabel__OBAtom_p_t__iterator(openbabelJNI.OBResidue_BeginAtoms(swigCPtr, this), true);
  }

  public SWIGTYPE_p_std__vectorTOpenBabel__OBAtom_p_t__iterator EndAtoms() {
    return new SWIGTYPE_p_std__vectorTOpenBabel__OBAtom_p_t__iterator(openbabelJNI.OBResidue_EndAtoms(swigCPtr, this), true);
  }

  public OBAtom BeginAtom(SWIGTYPE_p_std__vectorTOpenBabel__OBAtom_p_t__iterator i) {
    long cPtr = openbabelJNI.OBResidue_BeginAtom(swigCPtr, this, SWIGTYPE_p_std__vectorTOpenBabel__OBAtom_p_t__iterator.getCPtr(i));
    return (cPtr == 0) ? null : new OBAtom(cPtr, false);
  }

  public OBAtom NextAtom(SWIGTYPE_p_std__vectorTOpenBabel__OBAtom_p_t__iterator i) {
    long cPtr = openbabelJNI.OBResidue_NextAtom(swigCPtr, this, SWIGTYPE_p_std__vectorTOpenBabel__OBAtom_p_t__iterator.getCPtr(i));
    return (cPtr == 0) ? null : new OBAtom(cPtr, false);
  }

}
