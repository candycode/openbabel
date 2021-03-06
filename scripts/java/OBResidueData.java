/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 1.3.31
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */


public class OBResidueData extends OBGlobalDataBase {
  private long swigCPtr;

  protected OBResidueData(long cPtr, boolean cMemoryOwn) {
    super(openbabelJNI.SWIGOBResidueDataUpcast(cPtr), cMemoryOwn);
    swigCPtr = cPtr;
  }

  protected static long getCPtr(OBResidueData obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if(swigCPtr != 0 && swigCMemOwn) {
      swigCMemOwn = false;
      openbabelJNI.delete_OBResidueData(swigCPtr);
    }
    swigCPtr = 0;
    super.delete();
  }

  public OBResidueData() {
    this(openbabelJNI.new_OBResidueData(), true);
  }

  public void ParseLine(String arg0) {
    openbabelJNI.OBResidueData_ParseLine(swigCPtr, this, arg0);
  }

  public long GetSize() {
    return openbabelJNI.OBResidueData_GetSize(swigCPtr, this);
  }

  public boolean SetResName(String arg0) {
    return openbabelJNI.OBResidueData_SetResName(swigCPtr, this, arg0);
  }

  public int LookupBO(String arg0) {
    return openbabelJNI.OBResidueData_LookupBO__SWIG_0(swigCPtr, this, arg0);
  }

  public int LookupBO(String arg0, String arg1) {
    return openbabelJNI.OBResidueData_LookupBO__SWIG_1(swigCPtr, this, arg0, arg1);
  }

  public boolean LookupType(String arg0, SWIGTYPE_p_std__string arg1, SWIGTYPE_p_int arg2) {
    return openbabelJNI.OBResidueData_LookupType(swigCPtr, this, arg0, SWIGTYPE_p_std__string.getCPtr(arg1), SWIGTYPE_p_int.getCPtr(arg2));
  }

  public boolean AssignBonds(OBMol arg0, OBBitVec arg1) {
    return openbabelJNI.OBResidueData_AssignBonds(swigCPtr, this, OBMol.getCPtr(arg0), arg0, OBBitVec.getCPtr(arg1), arg1);
  }

}
