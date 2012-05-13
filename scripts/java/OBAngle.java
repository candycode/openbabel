/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 1.3.31
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */


public class OBAngle {
  private long swigCPtr;
  protected boolean swigCMemOwn;

  protected OBAngle(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(OBAngle obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if(swigCPtr != 0 && swigCMemOwn) {
      swigCMemOwn = false;
      openbabelJNI.delete_OBAngle(swigCPtr);
    }
    swigCPtr = 0;
  }

  public OBAngle(OBAngle arg0) {
    this(openbabelJNI.new_OBAngle__SWIG_2(OBAngle.getCPtr(arg0), arg0), true);
  }

  public void Clear() {
    openbabelJNI.OBAngle_Clear(swigCPtr, this);
  }

  public double GetAngle() {
    return openbabelJNI.OBAngle_GetAngle(swigCPtr, this);
  }

  public void SetAngle(double angle) {
    openbabelJNI.OBAngle_SetAngle(swigCPtr, this, angle);
  }

  public void SetAtoms(OBAtom vertex, OBAtom a, OBAtom b) {
    openbabelJNI.OBAngle_SetAtoms__SWIG_0(swigCPtr, this, OBAtom.getCPtr(vertex), vertex, OBAtom.getCPtr(a), a, OBAtom.getCPtr(b), b);
  }

  public void SetAtoms(SWIGTYPE_p_OpenBabel__tripleTOpenBabel__OBAtom_p_OpenBabel__OBAtom_p_OpenBabel__OBAtom_p_t atoms) {
    openbabelJNI.OBAngle_SetAtoms__SWIG_1(swigCPtr, this, SWIGTYPE_p_OpenBabel__tripleTOpenBabel__OBAtom_p_OpenBabel__OBAtom_p_OpenBabel__OBAtom_p_t.getCPtr(atoms));
  }

}
