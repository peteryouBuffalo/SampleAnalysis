// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME Processor_dict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "src/Processor.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_Processor(void *p = 0);
   static void *newArray_Processor(Long_t size, void *p);
   static void delete_Processor(void *p);
   static void deleteArray_Processor(void *p);
   static void destruct_Processor(void *p);
   static void streamer_Processor(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Processor*)
   {
      ::Processor *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::Processor >(0);
      static ::ROOT::TGenericClassInfo 
         instance("Processor", ::Processor::Class_Version(), "src/Processor.h", 22,
                  typeid(::Processor), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::Processor::Dictionary, isa_proxy, 16,
                  sizeof(::Processor) );
      instance.SetNew(&new_Processor);
      instance.SetNewArray(&newArray_Processor);
      instance.SetDelete(&delete_Processor);
      instance.SetDeleteArray(&deleteArray_Processor);
      instance.SetDestructor(&destruct_Processor);
      instance.SetStreamerFunc(&streamer_Processor);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Processor*)
   {
      return GenerateInitInstanceLocal((::Processor*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::Processor*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr Processor::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Processor::Class_Name()
{
   return "Processor";
}

//______________________________________________________________________________
const char *Processor::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Processor*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Processor::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Processor*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Processor::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Processor*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Processor::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Processor*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void Processor::Streamer(TBuffer &R__b)
{
   // Stream an object of class Processor.

   Reader::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_Processor(void *p) {
      return  p ? new(p) ::Processor : new ::Processor;
   }
   static void *newArray_Processor(Long_t nElements, void *p) {
      return p ? new(p) ::Processor[nElements] : new ::Processor[nElements];
   }
   // Wrapper around operator delete
   static void delete_Processor(void *p) {
      delete ((::Processor*)p);
   }
   static void deleteArray_Processor(void *p) {
      delete [] ((::Processor*)p);
   }
   static void destruct_Processor(void *p) {
      typedef ::Processor current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_Processor(TBuffer &buf, void *obj) {
      ((::Processor*)obj)->::Processor::Streamer(buf);
   }
} // end of namespace ROOT for class ::Processor

namespace {
  void TriggerDictionaryInitialization_Processor_dict_Impl() {
    static const char* headers[] = {
"src/Processor.h",
0
    };
    static const char* includePaths[] = {
"/cvmfs/cms.cern.ch/slc7_amd64_gcc820/lcg/root/6.14.09/include",
"/uscms_data/d3/peteryou/CMSSW_10_6_4/src/VHanalysis/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "Processor_dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$src/Processor.h")))  Processor;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "Processor_dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "src/Processor.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"Processor", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("Processor_dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_Processor_dict_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_Processor_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_Processor_dict() {
  TriggerDictionaryInitialization_Processor_dict_Impl();
}
