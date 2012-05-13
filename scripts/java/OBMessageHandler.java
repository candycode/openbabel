/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 1.3.31
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */


public class OBMessageHandler {
  private long swigCPtr;
  protected boolean swigCMemOwn;

  protected OBMessageHandler(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(OBMessageHandler obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if(swigCPtr != 0 && swigCMemOwn) {
      swigCMemOwn = false;
      openbabelJNI.delete_OBMessageHandler(swigCPtr);
    }
    swigCPtr = 0;
  }

  public OBMessageHandler() {
    this(openbabelJNI.new_OBMessageHandler(), true);
  }

  public void ThrowError(OBError err) {
    openbabelJNI.OBMessageHandler_ThrowError__SWIG_0(swigCPtr, this, OBError.getCPtr(err), err);
  }

  public void ThrowError(String method, String errorMsg, obMessageLevel level) {
    openbabelJNI.OBMessageHandler_ThrowError__SWIG_1(swigCPtr, this, method, errorMsg, level.swigValue());
  }

  public void ThrowError(String method, String errorMsg) {
    openbabelJNI.OBMessageHandler_ThrowError__SWIG_2(swigCPtr, this, method, errorMsg);
  }

  public vectorString GetMessagesOfLevel(obMessageLevel arg0) {
    return new vectorString(openbabelJNI.OBMessageHandler_GetMessagesOfLevel(swigCPtr, this, arg0.swigValue()), true);
  }

  public void StartLogging() {
    openbabelJNI.OBMessageHandler_StartLogging(swigCPtr, this);
  }

  public void StopLogging() {
    openbabelJNI.OBMessageHandler_StopLogging(swigCPtr, this);
  }

  public void SetMaxLogEntries(long max) {
    openbabelJNI.OBMessageHandler_SetMaxLogEntries(swigCPtr, this, max);
  }

  public long GetMaxLogEntries() {
    return openbabelJNI.OBMessageHandler_GetMaxLogEntries(swigCPtr, this);
  }

  public void ClearLog() {
    openbabelJNI.OBMessageHandler_ClearLog(swigCPtr, this);
  }

  public void SetOutputLevel(obMessageLevel level) {
    openbabelJNI.OBMessageHandler_SetOutputLevel(swigCPtr, this, level.swigValue());
  }

  public obMessageLevel GetOutputLevel() {
    return obMessageLevel.swigToEnum(openbabelJNI.OBMessageHandler_GetOutputLevel(swigCPtr, this));
  }

  public void SetOutputStream(SWIGTYPE_p_std__ostream os) {
    openbabelJNI.OBMessageHandler_SetOutputStream(swigCPtr, this, SWIGTYPE_p_std__ostream.getCPtr(os));
  }

  public SWIGTYPE_p_std__ostream GetOutputStream() {
    long cPtr = openbabelJNI.OBMessageHandler_GetOutputStream(swigCPtr, this);
    return (cPtr == 0) ? null : new SWIGTYPE_p_std__ostream(cPtr, false);
  }

  public boolean StartErrorWrap() {
    return openbabelJNI.OBMessageHandler_StartErrorWrap(swigCPtr, this);
  }

  public boolean StopErrorWrap() {
    return openbabelJNI.OBMessageHandler_StopErrorWrap(swigCPtr, this);
  }

  public long GetErrorMessageCount() {
    return openbabelJNI.OBMessageHandler_GetErrorMessageCount(swigCPtr, this);
  }

  public long GetWarningMessageCount() {
    return openbabelJNI.OBMessageHandler_GetWarningMessageCount(swigCPtr, this);
  }

  public long GetInfoMessageCount() {
    return openbabelJNI.OBMessageHandler_GetInfoMessageCount(swigCPtr, this);
  }

  public long GetAuditMessageCount() {
    return openbabelJNI.OBMessageHandler_GetAuditMessageCount(swigCPtr, this);
  }

  public long GetDebugMessageCount() {
    return openbabelJNI.OBMessageHandler_GetDebugMessageCount(swigCPtr, this);
  }

  public String GetMessageSummary() {
    return openbabelJNI.OBMessageHandler_GetMessageSummary(swigCPtr, this);
  }

}
