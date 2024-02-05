// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME Reader_dict

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
#include "src/Reader.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_Reader(void *p = 0);
   static void *newArray_Reader(Long_t size, void *p);
   static void delete_Reader(void *p);
   static void deleteArray_Reader(void *p);
   static void destruct_Reader(void *p);
   static void streamer_Reader(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Reader*)
   {
      ::Reader *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::Reader >(0);
      static ::ROOT::TGenericClassInfo 
         instance("Reader", ::Reader::Class_Version(), "src/Reader.h", 22,
                  typeid(::Reader), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::Reader::Dictionary, isa_proxy, 16,
                  sizeof(::Reader) );
      instance.SetNew(&new_Reader);
      instance.SetNewArray(&newArray_Reader);
      instance.SetDelete(&delete_Reader);
      instance.SetDeleteArray(&deleteArray_Reader);
      instance.SetDestructor(&destruct_Reader);
      instance.SetStreamerFunc(&streamer_Reader);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Reader*)
   {
      return GenerateInitInstanceLocal((::Reader*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::Reader*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr Reader::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Reader::Class_Name()
{
   return "Reader";
}

//______________________________________________________________________________
const char *Reader::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Reader*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Reader::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Reader*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Reader::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Reader*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Reader::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Reader*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void Reader::Streamer(TBuffer &R__b)
{
   // Stream an object of class Reader.

   TSelector::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_Reader(void *p) {
      return  p ? new(p) ::Reader : new ::Reader;
   }
   static void *newArray_Reader(Long_t nElements, void *p) {
      return p ? new(p) ::Reader[nElements] : new ::Reader[nElements];
   }
   // Wrapper around operator delete
   static void delete_Reader(void *p) {
      delete ((::Reader*)p);
   }
   static void deleteArray_Reader(void *p) {
      delete [] ((::Reader*)p);
   }
   static void destruct_Reader(void *p) {
      typedef ::Reader current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_Reader(TBuffer &buf, void *obj) {
      ((::Reader*)obj)->::Reader::Streamer(buf);
   }
} // end of namespace ROOT for class ::Reader

namespace {
  void TriggerDictionaryInitialization_Reader_dict_Impl() {
    static const char* headers[] = {
"src/Reader.h",
0
    };
    static const char* includePaths[] = {
"/cvmfs/cms.cern.ch/slc7_amd64_gcc820/lcg/root/6.14.09/include",
"/uscms_data/d3/peteryou/CMSSW_10_6_4/src/VHanalysis/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "Reader_dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$src/Reader.h")))  Reader;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "Reader_dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "src/Reader.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"Reader", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("Reader_dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_Reader_dict_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_Reader_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_Reader_dict() {
  TriggerDictionaryInitialization_Reader_dict_Impl();
}
