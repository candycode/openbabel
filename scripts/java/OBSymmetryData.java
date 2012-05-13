/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 1.3.31
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */


public class OBSymmetryData extends OBGenericData {
  private long swigCPtr;

  protected OBSymmetryData(long cPtr, boolean cMemoryOwn) {
    super(openbabelJNI.SWIGOBSymmetryDataUpcast(cPtr), cMemoryOwn);
    swigCPtr = cPtr;
  }

  protected static long getCPtr(OBSymmetryData obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if(swigCPtr != 0 && swigCMemOwn) {
      swigCMemOwn = false;
      openbabelJNI.delete_OBSymmetryData(swigCPtr);
    }
    swigCPtr = 0;
    super.delete();
  }

  public OBSymmetryData() {
    this(openbabelJNI.new_OBSymmetryData__SWIG_0(), true);
  }

  public OBSymmetryData(OBSymmetryData arg0) {
    this(openbabelJNI.new_OBSymmetryData__SWIG_1(OBSymmetryData.getCPtr(arg0), arg0), true);
  }

  public OBGenericData Clone(OBBase arg0) {
    long cPtr = openbabelJNI.OBSymmetryData_Clone(swigCPtr, this, OBBase.getCPtr(arg0), arg0);
    return (cPtr == 0) ? null : new OBGenericData(cPtr, false);
  }

  public void SetData(String pg, String sg) {
    openbabelJNI.OBSymmetryData_SetData__SWIG_0(swigCPtr, this, pg, sg);
  }

  public void SetData(String pg) {
    openbabelJNI.OBSymmetryData_SetData__SWIG_1(swigCPtr, this, pg);
  }

  public void SetPointGroup(String pg) {
    openbabelJNI.OBSymmetryData_SetPointGroup(swigCPtr, this, pg);
  }

  public void SetSpaceGroup(String sg) {
    openbabelJNI.OBSymmetryData_SetSpaceGroup(swigCPtr, this, sg);
  }

  public String GetPointGroup() {
    return openbabelJNI.OBSymmetryData_GetPointGroup(swigCPtr, this);
  }

  public String GetSpaceGroup() {
    return openbabelJNI.OBSymmetryData_GetSpaceGroup(swigCPtr, this);
  }

}
