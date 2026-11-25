////////////////////////////////////////////////////////////////////////
// SecretStrings.h
// (c) 2026 Wei-Min Chen
//
// This file is distributed as part of the SECRET source code package
// and may not be redistributed in any form, without prior written
// permission from the author. Permission is granted for you to
// modify this file for your own personal use, but modified versions
// must retain this copyright notice and must not be distributed.
//
// Permission is granted for you to use this file to compile SECRET.
//
// All computer programs have bugs. Use this file at your own risk.
//
// The text classes of SECRET, written for SECRET:
//
//   String          a text string with the usual operators; converts to const char *,
//                   int and double.  Comparisons order digit runs numerically
//                   ("x9" < "x10"); String::caseSensitive = false makes Compare ignore case.
//   StringArray     a list of strings; AddTokens / AddColumns split a line.  An element
//                   keeps its address while other elements are added.
//   StringIntHash   text -> int, keys are compared without regard to case.  The index
//                   that Add / Find return stays valid until that key is deleted.
//   StringMap       text -> pointer, kept in order; FindStem finds a unique abbreviation
//   StringIntMap    text -> int, the same
//
// The code is C++98 so that older compilers (C++ Builder 6) can build it.
////////////////////////////////////////////////////////////////////////
#ifndef __SECRETSTRINGS_H__
#define __SECRETSTRINGS_H__

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string>
#include <vector>
#include <map>

class String
   {
   protected:
      std::string text;

   public:
      static bool caseSensitive;       // used by Compare, Find and FindChar

      String(int startsize = 0)        { if (startsize > 0) text.reserve(startsize); }
      String(const char * s)           { if (s != NULL) text = s; }
      String(const String & s)         : text(s.text) { }
      String(char ch, int count = 1)   { if (count > 0) text.assign(count, ch); }

      String & Clear()                 { text.erase(); return *this; }
      bool     IsEmpty() const         { return text.empty(); }
      int      Length() const          { return (int) text.size(); }

      String & Copy(const String & s)  { text = s.text; return *this; }
      String & Copy(const String & s, int start, int count);
      String & Copy(const char * s)    { text = (s == NULL) ? "" : s; return *this; }

      void Swap(String & s)            { text.swap(s.text); }

      String & ToUpper();
      String & ToLower();
      String   AsUpper() const         { String t(*this); return t.ToUpper(); }
      String   AsLower() const         { String t(*this); return t.ToLower(); }

      // assignment and concatenation; numbers are written as text: int and unsigned
      // int in decimal, double with three decimals
      String & operator =  (char ch);
      String   operator +  (char ch) const         { String r(*this); r += ch; return r; }
      String & operator += (char ch);

      String & operator =  (const String & rhs)    { text = rhs.text; return *this; }
      String   operator +  (const String & rhs) const { String r(*this); r.text += rhs.text; return r; }
      String & operator += (const String & rhs)    { text += rhs.text; return *this; }

      String & operator =  (const char * rhs)      { return Copy(rhs); }
      String   operator +  (const char * rhs) const { String r(*this); r += rhs; return r; }
      String & operator += (const char * rhs)      { if (rhs != NULL) text += rhs; return *this; }

      String & operator =  (int rhs);
      String   operator +  (int rhs) const         { String r(*this); r += rhs; return r; }
      String & operator += (int rhs);

      String & operator =  (unsigned int rhs);
      String   operator +  (unsigned int rhs) const { String r(*this); r += rhs; return r; }
      String & operator += (unsigned int rhs);

      String & operator =  (double rhs);
      String   operator +  (double rhs) const      { String r(*this); r += rhs; return r; }
      String & operator += (double rhs);

      String   operator *  (unsigned int rhs) const { String r(*this); r *= rhs; return r; }
      String & operator *= (unsigned int rhs);     // repeat the text rhs times

      // < 0, 0, > 0 like strcmp, but a run of digits is compared as a number when the
      // strings first differ inside it
      int Compare(const String & rhs) const        { return caseSensitive ? FastCompare(rhs) : SlowCompare(rhs); }
      int Compare(const char * rhs) const          { return caseSensitive ? FastCompare(rhs) : SlowCompare(rhs); }
      int FastCompare(const String & rhs) const    { return CompareText(text.c_str(), rhs.text.c_str(), false); }
      int FastCompare(const char * rhs) const;
      int SlowCompare(const String & rhs) const    { return CompareText(text.c_str(), rhs.text.c_str(), true); }
      int SlowCompare(const char * rhs) const;

      // compare the beginning of this string with the stem (0: this starts with the stem)
      int CompareToStem(const String & stem) const { return caseSensitive ? FastCompareToStem(stem) : SlowCompareToStem(stem); }
      int FastCompareToStem(const String & stem) const { return CompareStem(text.c_str(), stem.text.c_str(), false); }
      int SlowCompareToStem(const String & stem) const { return CompareStem(text.c_str(), stem.text.c_str(), true); }

      int operator == (const String & rhs) const   { return Compare(rhs) == 0; }
      int operator != (const String & rhs) const   { return Compare(rhs) != 0; }
      int operator <  (const String & rhs) const   { return Compare(rhs)  < 0; }
      int operator >  (const String & rhs) const   { return Compare(rhs)  > 0; }
      int operator >= (const String & rhs) const   { return Compare(rhs) >= 0; }
      int operator <= (const String & rhs) const   { return Compare(rhs) <= 0; }

      int operator == (const char * rhs) const     { return Compare(rhs) == 0; }
      int operator != (const char * rhs) const     { return Compare(rhs) != 0; }
      int operator <  (const char * rhs) const     { return Compare(rhs)  < 0; }
      int operator >  (const char * rhs) const     { return Compare(rhs)  > 0; }
      int operator <= (const char * rhs) const     { return Compare(rhs) <= 0; }
      int operator >= (const char * rhs) const     { return Compare(rhs) >= 0; }

      operator const char * () const   { return text.c_str(); }
      operator int () const            { return atoi(text.c_str()); }
      operator double () const         { return atof(text.c_str()); }
      char   operator [] (int i) const { return text.c_str()[i]; }     // [Length()] is the end mark
      char & operator [] (int i)       { return text[i]; }

      char & Last()                    { return text[text.size() - 1]; }
      char & First()                   { return text[0]; }

      String Left(int count) const;
      String Right(int count) const;
      String Mid(int start, int end) const;             // characters start..end (both included)
      String SubStr(int start, int count) const;
      String SubStr(int start) const;

      int FindChar(char ch, int start = 0) const       { return caseSensitive ? FastFindChar(ch, start) : SlowFindChar(ch, start); }
      int FastFindChar(char ch, int start = 0) const;
      int SlowFindChar(char ch, int start = 0) const;

      int Find(const String & str, int start = 0) const  { return caseSensitive ? FastFind(str, start) : SlowFind(str, start); }
      int FastFind(const String & str, int start = 0) const;
      int SlowFind(const String & str, int start = 0) const;

      String & ExcludeCharacters(const String & s);     // remove every character that is in s
      String & ExcludeCharacters(const char * s)        { return ExcludeCharacters(String(s)); }

      int SetLength(int newlen);                        // cut, or pad with blanks
      int Dimension(int newlen)        { return SetLength(newlen); }

      String & Add(const String & s)   { return *this += s; }
      String & Add(char ch)            { return *this += ch; }

      String & Invert();                                // reverse the text
      String & Trim();                                  // remove white space at both ends
      String & Trim(char character);                    // remove one character at both ends

      long   AsInteger() const;                         // decimal, or hexadecimal with 0x
      double AsDouble() const          { return atof(text.c_str()); }
      bool   IsNumber() const;

      int printf(const char * format, ...);
      int vprintf(const char * format, va_list arglist);
      int catprintf(const char * format, ...);
      int vcatprintf(const char * format, va_list arglist);

      const unsigned char * uchar() const { return (const unsigned char *) text.c_str(); }

   private:
      static int CompareText(const char * a, const char * b, bool ignoreCase);
      static int CompareStem(const char * a, const char * stem, bool ignoreCase);
   };

inline int Compare(const String & s1, const String & s2)     { return s1.Compare(s2); }
inline int Compare(const String & s1, const char * s2)       { return s1.Compare(s2); }
inline int Compare(const char * s1, const String & s2)       { return -s2.Compare(s1); }

String operator + (char lhs, const String & rhs);
String operator + (const char * lhs, const String & rhs);
String operator + (int lhs, const String & rhs);
String operator + (unsigned int lhs, const String & rhs);

class StringArray
   {
   protected:
      std::vector<String *> items;       // pointers: an element keeps its address while others are added

   public:
      StringArray(int startsize = 0);
      StringArray(const StringArray & original);
      ~StringArray();

      // text files: one element per line
      void Write(FILE * f) const;
      void Write(const char * filename) const;
      void WriteLine(FILE * f) const;                    // all elements on one line, tab-delimited
      void Print() const;
      void PrintLine() const;

      void Clear();

      int Length() const               { return (int) items.size(); }
      int Dimension(int newcount);                       // new elements are empty strings
      int CharLength() const;

      String & operator [] (int i)             { return *(items[i]); }
      const String & operator [] (int i) const { return *(items[i]); }

      // split a string and append the pieces; the result is the new length of the array
      int AddColumns(const String & s, char ch = '\t');  // every separator counts: empty pieces are kept
      int AddTokens(const String & s, char ch);          // runs of the separator count once
      int AddTokens(const String & s, const String & separators = " \t\r\n");

      int ReplaceColumns(const String & s, char ch = '\t')
         { Clear(); return AddColumns(s, ch); }
      int ReplaceTokens(const String & s, const String & separators = " \t\r\n")
         { Clear(); return AddTokens(s, separators); }

      int  Add(const String & s);                        // the result is the new length
      void InsertAt(int position, const String & s);
      void Delete(int position);

      String & Last() const;                             // stops with an error when empty
      int      Push(const String & s)  { return Add(s); }
      String   Pop();

      // linear search; -1 when absent
      int Find(const String & s) const;
      int FastFind(const String & s) const;
      int SlowFind(const String & s) const;

      void Sort();                                       // with Compare
      void Trim();                                       // Trim() of every element

      StringArray & operator = (const StringArray & rhs);
      bool operator == (const StringArray & rhs) const;
      bool operator != (const StringArray & rhs) const  { return !(*this == rhs); }
   };

// text -> int; the case of a key does not matter
class StringIntHash
   {
   protected:
      std::vector<String *> keys;        // a slot: key (NULL after Delete) and value
      std::vector<int>      values;
      std::map<std::string, int> slotOf; // upper case key -> slot
      int count;

      static std::string Fold(const String & key);

   public:
      StringIntHash(int startsize = 32);
      ~StringIntHash();

      void Clear();

      int  Capacity() const            { return (int) keys.size(); }   // slots to look at with SlotInUse
      int  Entries() const             { return count; }
      bool SlotInUse(int index) const  { return keys[index] != NULL; }

      int Integer(int i) const         { return values[i]; }
      int Integer(const String & key) const;                           // -1 when absent
      void SetInteger(int i, int value)               { values[i] = value; }
      void SetInteger(const String & key, int value)  { Add(key, value); }

      int IncrementCount(const String & key);          // absent counts as 0
      int DecrementCount(const String & key);
      int GetCount(const String & key) const;          // 0 when absent
      int GetCount(int index) const    { return values[index]; }

      int Add(const String & s, int integer);          // adds or changes; the result is the index
      int Find(const String & s, int defaultValue);    // adds with the default when absent
      int Find(const String & s) const;                // -1 when absent

      StringIntHash & operator = (const StringIntHash & rhs);

      const String & operator [] (int i) const        { return *(keys[i]); }
      String & operator [] (int i)                    { return *(keys[i]); }

      void Delete(int index);
      void Delete(const String & key)                 { Delete(Find(key)); }

   private:
      StringIntHash(const StringIntHash &);            // not copyable: use operator =
   };

// text -> pointer, in order (case does not matter).  The index of a key changes when
// another key is added before it.
class StringMap
   {
   protected:
      std::vector<String *> keys;
      std::vector<void *>   objects;

   public:
      StringMap(int startsize = 0);
      ~StringMap();

      void Clear();
      int  Length() const              { return (int) keys.size(); }

      void * Object(int i) const       { return objects[i]; }
      void * Object(const String & key) const
         {
         int index = Find(key);
         return (index >= 0) ? objects[index] : NULL;
         }
      void SetObject(int i, void * object)             { objects[i] = object; }
      void SetObject(const String & key, void * object) { Add(key, object); }

      int Add(const String & s, void * object = NULL); // adds or changes; the result is the index
      int Find(const String & s) const;                // -1 when absent
      // the key that starts with the stem: -1 for none, -2 for more than one
      int FindStem(const String & stem) const;

      const String & operator [] (int i) const        { return *(keys[i]); }
      String & operator [] (int i)                    { return *(keys[i]); }

      void Delete(int index);

   private:
      StringMap(const StringMap &);
      StringMap & operator = (const StringMap &);
      int Locate(const String & s, bool & found) const;
   };

// text -> int, in order (case does not matter)
class StringIntMap
   {
   protected:
      std::vector<String *> keys;
      std::vector<int>      values;

   public:
      StringIntMap(int startsize = 0);
      ~StringIntMap();

      void Clear();
      int  Length() const              { return (int) keys.size(); }

      int Integer(int i) const         { return values[i]; }
      int Integer(const String & key) const;           // -1 when absent
      void SetInteger(int i, int value)               { values[i] = value; }
      void SetInteger(const String & key, int value)  { Add(key, value); }

      int Add(const String & s, int i);
      int Find(const String & s, int defaultValue);    // adds with the default when absent
      int Find(const String & s) const;
      int FindStem(const String & stem) const;

      int IncrementCount(const String & key);          // absent counts as 0
      int DecrementCount(const String & key);
      int GetCount(const String & key) const;          // 0 when absent
      int GetCount(int index) const    { return values[index]; }

      const String & operator [] (int i) const        { return *(keys[i]); }
      String & operator [] (int i)                    { return *(keys[i]); }

      void Delete(int index);

   private:
      StringIntMap(const StringIntMap &);
      StringIntMap & operator = (const StringIntMap &);
      int Locate(const String & s, bool & found) const;
   };

#endif
