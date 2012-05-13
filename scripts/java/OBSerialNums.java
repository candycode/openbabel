/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 1.3.31
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */


public class OBSerialNums extends OBGenericData {
  private long swigCPtr;

  protected OBSerialNums(long cPtr, boolean cMemoryOwn) {
    super(openbabelJNI.SWIGOBSerialNumsUpcast(cPtr), cMemoryOwn);
    swigCPtr = cPtr;
  }

  protected static long getCPtr(OBSerialNums obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if(swigCPtr != 0 && swigCMemOwn) {
      swigCMemOwn = false;
      openbabelJNI.delete_OBSerialNums(swigCPtr);
    }
    swigCPtr = 0;
    super.delete();
  }

  public OBSerialNums() {
    this(openbabelJNI.new_OBSerialNums__SWIG_0(), true);
  }

  public OBSerialNums(OBSerialNums cp) {
    this(openbabelJNI.new_OBSerialNums__SWIG_1(OBSerialNums.getCPtr(cp), cp), true);
  }

  public OBGenericData Clone(OBBase arg0) {
    long cPtr = openbabelJNI.OBSerialNums_Clone(swigCPtr, this, OBBase.getCPtr(arg0), arg0);
    return (cPtr == 0) ? null : new OBGenericData(cPtr, false);
  }

  public SWIGTYPE_p_std__mapTint_OpenBabel__OBAtom_p_t GetData() {
    return new SWIGTYPE_p_std__mapTint_OpenBabel__OBAtom_p_t(openbabelJNI.OBSerialNums_GetData(swigCPtr, this), false);
  }

  public void SetData(SWIGTYPE_p_std__mapTint_OpenBabel__OBAtom_p_t sm) {
    openbabelJNI.OBSerialNums_SetData(swigCPtr, this, SWIGTYPE_p_std__mapTint_OpenBabel__OBAtom_p_t.getCPtr(sm));
  }

}
