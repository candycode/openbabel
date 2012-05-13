/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 1.3.31
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */


public class OBTypeTable extends OBGlobalDataBase {
  private long swigCPtr;

  protected OBTypeTable(long cPtr, boolean cMemoryOwn) {
    super(openbabelJNI.SWIGOBTypeTableUpcast(cPtr), cMemoryOwn);
    swigCPtr = cPtr;
  }

  protected static long getCPtr(OBTypeTable obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if(swigCPtr != 0 && swigCMemOwn) {
      swigCMemOwn = false;
      openbabelJNI.delete_OBTypeTable(swigCPtr);
    }
    swigCPtr = 0;
    super.delete();
  }

  public OBTypeTable() {
    this(openbabelJNI.new_OBTypeTable(), true);
  }

  public void ParseLine(String arg0) {
    openbabelJNI.OBTypeTable_ParseLine(swigCPtr, this, arg0);
  }

  public long GetSize() {
    return openbabelJNI.OBTypeTable_GetSize(swigCPtr, this);
  }

  public boolean SetFromType(String arg0) {
    return openbabelJNI.OBTypeTable_SetFromType(swigCPtr, this, arg0);
  }

  public boolean SetToType(String arg0) {
    return openbabelJNI.OBTypeTable_SetToType(swigCPtr, this, arg0);
  }

  public boolean Translate(String to, String from) {
    return openbabelJNI.OBTypeTable_Translate__SWIG_0(swigCPtr, this, to, from);
  }

  public boolean Translate(SWIGTYPE_p_std__string to, String from) {
    return openbabelJNI.OBTypeTable_Translate__SWIG_1(swigCPtr, this, SWIGTYPE_p_std__string.getCPtr(to), from);
  }

  public String Translate(String from) {
    return openbabelJNI.OBTypeTable_Translate__SWIG_2(swigCPtr, this, from);
  }

  public String GetFromType() {
    return openbabelJNI.OBTypeTable_GetFromType(swigCPtr, this);
  }

  public String GetToType() {
    return openbabelJNI.OBTypeTable_GetToType(swigCPtr, this);
  }

}
