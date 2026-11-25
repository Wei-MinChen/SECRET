////////////////////////////////////////////////////////////////////////
// SecretStrings.cpp
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
// The text classes of SECRET; see SecretStrings.h.
////////////////////////////////////////////////////////////////////////

#include "SecretStrings.h"
#include "Error.h"

#include <ctype.h>
#include <string.h>
#include <algorithm>

bool String::caseSensitive = true;

//////////////////////////////////////////////////////////////////////
// String
//////////////////////////////////////////////////////////////////////

// Compare two C strings.  At the first difference the strings are ordered by that
// character, unless it is inside a run of digits: then the string with the longer
// run of digits is the larger (so "x9" < "x10").
int String::CompareText(const char * a, const char * b, bool ignoreCase)
   {
   for (int i = 0; ; i++)
      {
      int x = a[i], y = b[i];

      if (ignoreCase)
         {
         x = toupper(x);
         y = toupper(y);
         }

      if (x != y)
         {
         int d = i;

         while (isdigit((unsigned char) a[d]) && isdigit((unsigned char) b[d]))
            d++;

         if (isdigit((unsigned char) a[d])) return 1;
         if (isdigit((unsigned char) b[d])) return -1;

         return x - y;
         }

      if (a[i] == 0) return 0;
      }
   }

// 0 when a starts with the stem
int String::CompareStem(const char * a, const char * stem, bool ignoreCase)
   {
   for (int i = 0; stem[i] != 0; i++)
      {
      int x = a[i], y = stem[i];

      if (ignoreCase)
         {
         x = toupper(x);
         y = toupper(y);
         }

      if (x != y) return x - y;
      }

   return 0;
   }

int String::FastCompare(const char * rhs) const
   {
   if (rhs == NULL) return -Length();
   return CompareText(text.c_str(), rhs, false);
   }

int String::SlowCompare(const char * rhs) const
   {
   if (rhs == NULL) return -Length();
   return CompareText(text.c_str(), rhs, true);
   }

String & String::Copy(const String & s, int start, int count)
   {
   String piece = s.SubStr(start, count);
   text = piece.text;
   return *this;
   }

String & String::ToUpper()
   {
   for (size_t i = 0; i < text.size(); i++)
      text[i] = (char) toupper((unsigned char) text[i]);
   return *this;
   }

String & String::ToLower()
   {
   for (size_t i = 0; i < text.size(); i++)
      text[i] = (char) tolower((unsigned char) text[i]);
   return *this;
   }

String & String::operator = (char ch)
   {
   text.erase();
   if (ch) text.assign(1, ch);
   return *this;
   }

String & String::operator += (char ch)
   {
   if (ch) text.append(1, ch);
   return *this;
   }

String & String::operator = (int rhs)
   {
   char number[32];
   sprintf(number, "%d", rhs);
   text = number;
   return *this;
   }

String & String::operator += (int rhs)
   {
   char number[32];
   sprintf(number, "%d", rhs);
   text += number;
   return *this;
   }

String & String::operator = (unsigned int rhs)
   {
   char number[32];
   sprintf(number, "%u", rhs);
   text = number;
   return *this;
   }

String & String::operator += (unsigned int rhs)
   {
   char number[32];
   sprintf(number, "%u", rhs);
   text += number;
   return *this;
   }

String & String::operator = (double rhs)
   {
   char number[400];                     // "%.3f" of the largest double has 313 characters
   sprintf(number, "%.3f", rhs);
   text = number;
   return *this;
   }

String & String::operator += (double rhs)
   {
   String number;
   number = rhs;
   text += number.text;
   return *this;
   }

String & String::operator *= (unsigned int rhs)
   {
   if (rhs == 0)
      text.erase();
   else
      {
      std::string original(text);

      text.reserve(original.size() * rhs);
      for (unsigned int i = 1; i < rhs; i++)
         text += original;
      }
   return *this;
   }

String String::Left(int count) const
   {
   if (count < 0) count = 0;
   if (count > Length()) count = Length();

   String result;
   result.text = text.substr(0, count);
   return result;
   }

String String::Right(int count) const
   {
   if (count < 0) count = 0;
   if (count > Length()) count = Length();

   String result;
   result.text = text.substr(text.size() - count);
   return result;
   }

String String::SubStr(int start, int count) const
   {
   String result;

   if (start < 0)
      {
      count += start;
      start = 0;
      }

   if (start >= Length() || count <= 0)
      return result;

   if (count > Length() - start)
      count = Length() - start;

   result.text = text.substr(start, count);
   return result;
   }

String String::SubStr(int start) const
   {
   return SubStr(start, Length() - start);
   }

String String::Mid(int start, int end) const
   {
   return SubStr(start, end - start + 1);
   }

int String::FastFindChar(char ch, int start) const
   {
   if (start < 0) start = 0;

   for ( ; start < Length(); start++)
      if (text[start] == ch)
         return start;
   return -1;
   }

int String::SlowFindChar(char ch, int start) const
   {
   if (start < 0) start = 0;

   int upper = toupper((unsigned char) ch);

   for ( ; start < Length(); start++)
      if (toupper((unsigned char) text[start]) == upper)
         return start;
   return -1;
   }

// An empty pattern is never found
int String::FastFind(const String & str, int start) const
   {
   if (str.IsEmpty()) return -1;
   if (start < 0) start = 0;

   std::string::size_type at = text.find(str.text, start);

   return at == std::string::npos ? -1 : (int) at;
   }

int String::SlowFind(const String & str, int start) const
   {
   if (str.IsEmpty()) return -1;
   if (start < 0) start = 0;

   for ( ; start + str.Length() <= Length(); start++)
      {
      int i = 0;

      while (i < str.Length() && toupper((unsigned char) text[start + i]) ==
                                 toupper((unsigned char) str.text[i]))
         i++;

      if (i == str.Length()) return start;
      }
   return -1;
   }

String & String::ExcludeCharacters(const String & s)
   {
   std::string kept;

   for (size_t i = 0; i < text.size(); i++)
      if (s.FindChar(text[i]) == -1)
         kept += text[i];

   text = kept;
   return *this;
   }

int String::SetLength(int newlen)
   {
   if (newlen < 0) newlen = 0;
   text.resize(newlen, ' ');
   return newlen;
   }

String & String::Invert()
   {
   std::reverse(text.begin(), text.end());
   return *this;
   }

String & String::Trim()
   {
   size_t first = 0, last = text.size();

   while (first < last && isspace((unsigned char) text[first]))
      first++;
   while (last > first && isspace((unsigned char) text[last - 1]))
      last--;

   text = text.substr(first, last - first);
   return *this;
   }

String & String::Trim(char character)
   {
   size_t first = 0, last = text.size();

   while (first < last && text[first] == character)
      first++;
   while (last > first && text[last - 1] == character)
      last--;

   text = text.substr(first, last - first);
   return *this;
   }

long String::AsInteger() const
   {
   const char * p = text.c_str();
   long value = 0;
   int  base = 10, sign = 1;

   if (*p == '-')
      {
      sign = -1;
      p++;
      }

   if (p[0] == '0' && (p[1] == 'x' || p[1] == 'X') && p[2] != 0)
      {
      base = 16;
      p += 2;
      }

   for ( ; *p != 0; p++)
      {
      int digit = toupper((unsigned char) *p);

      if (digit >= '0' && digit <= '9')
         value = value * base + (digit - '0');
      else if (base == 16 && digit >= 'A' && digit <= 'F')
         value = value * base + (digit - 'A' + 10);
      else
         break;
      }

   return sign * value;
   }

// [+-] digits [. digits] [(e|E) [+-] digits], with at least one digit in the mantissa
bool String::IsNumber() const
   {
   const char * p = text.c_str();
   bool digits = false;

   if (*p == '-' || *p == '+') p++;

   while (*p >= '0' && *p <= '9') { p++; digits = true; }

   if (*p == '.')
      {
      p++;
      while (*p >= '0' && *p <= '9') { p++; digits = true; }
      }

   if (!digits) return false;

   if (*p == 'e' || *p == 'E')
      {
      p++;
      if (*p == '-' || *p == '+') p++;

      digits = false;
      while (*p >= '0' && *p <= '9') { p++; digits = true; }
      }

   return *p == 0 && digits;
   }

int String::printf(const char * format, ...)
   {
   va_list ap;
   va_start(ap, format);
   int length = vprintf(format, ap);
   va_end(ap);
   return length;
   }

int String::catprintf(const char * format, ...)
   {
   va_list ap;
   va_start(ap, format);
   int length = vcatprintf(format, ap);
   va_end(ap);
   return length;
   }

// Format into text, with a buffer that is made larger until the result fits
int String::vprintf(const char * format, va_list ap)
   {
   std::vector<char> buffer(256);

   for (;;)
      {
#ifdef va_copy
      va_list again;
      va_copy(again, ap);
      int needed = vsnprintf(&buffer[0], buffer.size(), format, again);
      va_end(again);
#else
      int needed = vsnprintf(&buffer[0], buffer.size(), format, ap);
#endif

      if (needed >= 0 && needed < (int) buffer.size())
         {
         text.assign(&buffer[0], needed);
         return needed;
         }

      buffer.resize(needed >= 0 ? needed + 1 : buffer.size() * 2);
      }
   }

int String::vcatprintf(const char * format, va_list ap)
   {
   String added;

   added.vprintf(format, ap);
   text += added.text;
   return Length();
   }

String operator + (const char * lhs, const String & rhs)
   {
   String result(lhs);
   result += rhs;
   return result;
   }

String operator + (char lhs, const String & rhs)
   {
   String result;
   result += lhs;
   result += rhs;
   return result;
   }

String operator + (int lhs, const String & rhs)
   {
   String result;
   result = lhs;
   result += rhs;
   return result;
   }

String operator + (unsigned int lhs, const String & rhs)
   {
   String result;
   result = lhs;
   result += rhs;
   return result;
   }

//////////////////////////////////////////////////////////////////////
// StringArray
//////////////////////////////////////////////////////////////////////

StringArray::StringArray(int startsize)
   {
   Dimension(startsize);
   }

StringArray::StringArray(const StringArray & original)
   {
   for (int i = 0; i < original.Length(); i++)
      items.push_back(new String(original[i]));
   }

StringArray::~StringArray()
   {
   Clear();
   }

void StringArray::Clear()
   {
   for (size_t i = 0; i < items.size(); i++)
      delete items[i];
   items.clear();
   }

int StringArray::Dimension(int newcount)
   {
   if (newcount < 0) newcount = 0;

   while ((int) items.size() > newcount)
      {
      delete items.back();
      items.pop_back();
      }

   while ((int) items.size() < newcount)
      items.push_back(new String);

   return newcount;
   }

int StringArray::CharLength() const
   {
   int total = 0;

   for (size_t i = 0; i < items.size(); i++)
      total += items[i]->Length();
   return total;
   }

void StringArray::Write(FILE * f) const
   {
   for (size_t i = 0; i < items.size(); i++)
      fprintf(f, "%s\n", (const char *) *items[i]);
   }

void StringArray::Write(const char * filename) const
   {
   FILE * f = fopen(filename, "wt");

   if (f == NULL) return;
   Write(f);
   fclose(f);
   }

void StringArray::WriteLine(FILE * f) const
   {
   for (size_t i = 0; i < items.size(); i++)
      fprintf(f, "%s%c", (const char *) *items[i], i + 1 == items.size() ? '\n' : '\t');
   }

void StringArray::Print() const
   {
   Write(stdout);
   }

void StringArray::PrintLine() const
   {
   WriteLine(stdout);
   }

// Every separator ends a column, so "a\t\tb" has three columns
int StringArray::AddColumns(const String & s, char ch)
   {
   int start = 0;

   for (;;)
      {
      int end = s.FindChar(ch, start);

      if (end == -1)
         {
         items.push_back(new String(s.SubStr(start)));
         break;
         }

      items.push_back(new String(s.SubStr(start, end - start)));
      start = end + 1;
      }

   return Length();
   }

int StringArray::AddTokens(const String & s, char ch)
   {
   int pos = 0, len = s.Length();

   while (pos < len)
      {
      while (pos < len && s[pos] == ch) pos++;

      int start = pos;

      while (pos < len && s[pos] != ch) pos++;

      if (pos > start)
         items.push_back(new String(s.SubStr(start, pos - start)));
      }

   return Length();
   }

int StringArray::AddTokens(const String & s, const String & separators)
   {
   int pos = 0, len = s.Length();

   while (pos < len)
      {
      while (pos < len && separators.FindChar(s[pos]) != -1) pos++;

      int start = pos;

      while (pos < len && separators.FindChar(s[pos]) == -1) pos++;

      if (pos > start)
         items.push_back(new String(s.SubStr(start, pos - start)));
      }

   return Length();
   }

int StringArray::Add(const String & s)
   {
   items.push_back(new String(s));
   return Length();
   }

void StringArray::InsertAt(int position, const String & s)
   {
   items.insert(items.begin() + position, new String(s));
   }

void StringArray::Delete(int position)
   {
   delete items[position];
   items.erase(items.begin() + position);
   }

String & StringArray::Last() const
   {
   if (items.empty()) error("StringArray: Null String Access");
   return *(items.back());
   }

String StringArray::Pop()
   {
   if (items.empty()) error("StringArray: Null String Access");

   String result(*items.back());

   delete items.back();
   items.pop_back();
   return result;
   }

int StringArray::Find(const String & s) const
   {
   for (size_t i = 0; i < items.size(); i++)
      if (*items[i] == s)
         return (int) i;
   return -1;
   }

int StringArray::FastFind(const String & s) const
   {
   for (size_t i = 0; i < items.size(); i++)
      if (items[i]->FastCompare(s) == 0)
         return (int) i;
   return -1;
   }

int StringArray::SlowFind(const String & s) const
   {
   for (size_t i = 0; i < items.size(); i++)
      if (items[i]->SlowCompare(s) == 0)
         return (int) i;
   return -1;
   }

static bool StringPointerBefore(const String * a, const String * b)
   {
   return a->Compare(*b) < 0;
   }

void StringArray::Sort()
   {
   std::stable_sort(items.begin(), items.end(), StringPointerBefore);
   }

void StringArray::Trim()
   {
   for (size_t i = 0; i < items.size(); i++)
      items[i]->Trim();
   }

StringArray & StringArray::operator = (const StringArray & rhs)
   {
   if (this != &rhs)
      {
      Clear();
      for (int i = 0; i < rhs.Length(); i++)
         items.push_back(new String(rhs[i]));
      }
   return *this;
   }

bool StringArray::operator == (const StringArray & rhs) const
   {
   if (Length() != rhs.Length()) return false;

   for (int i = 0; i < Length(); i++)
      if (*items[i] != *rhs.items[i])
         return false;
   return true;
   }

//////////////////////////////////////////////////////////////////////
// StringIntHash
//////////////////////////////////////////////////////////////////////

StringIntHash::StringIntHash(int)
   {
   count = 0;
   }

StringIntHash::~StringIntHash()
   {
   Clear();
   }

std::string StringIntHash::Fold(const String & key)
   {
   std::string folded((const char *) key);

   for (size_t i = 0; i < folded.size(); i++)
      folded[i] = (char) toupper((unsigned char) folded[i]);
   return folded;
   }

void StringIntHash::Clear()
   {
   for (size_t i = 0; i < keys.size(); i++)
      delete keys[i];

   keys.clear();
   values.clear();
   slotOf.clear();
   count = 0;
   }

int StringIntHash::Find(const String & s) const
   {
   std::map<std::string, int>::const_iterator at = slotOf.find(Fold(s));

   return at == slotOf.end() ? -1 : at->second;
   }

int StringIntHash::Add(const String & s, int integer)
   {
   std::string folded = Fold(s);
   std::map<std::string, int>::iterator at = slotOf.find(folded);

   if (at != slotOf.end())
      {
      values[at->second] = integer;
      return at->second;
      }

   int slot = (int) keys.size();

   keys.push_back(new String(s));
   values.push_back(integer);
   slotOf[folded] = slot;
   count++;
   return slot;
   }

int StringIntHash::Find(const String & s, int defaultValue)
   {
   int slot = Find(s);

   return slot >= 0 ? slot : Add(s, defaultValue);
   }

int StringIntHash::Integer(const String & key) const
   {
   int slot = Find(key);

   return slot >= 0 ? values[slot] : -1;
   }

int StringIntHash::GetCount(const String & key) const
   {
   int slot = Find(key);

   return slot >= 0 ? values[slot] : 0;
   }

int StringIntHash::IncrementCount(const String & key)
   {
   int slot = Find(key);

   if (slot < 0)
      {
      Add(key, 1);
      return 1;
      }

   return ++values[slot];
   }

int StringIntHash::DecrementCount(const String & key)
   {
   int slot = Find(key);

   if (slot < 0)
      {
      Add(key, -1);
      return -1;
      }

   return --values[slot];
   }

void StringIntHash::Delete(int index)
   {
   if (index < 0 || index >= (int) keys.size() || keys[index] == NULL)
      return;

   slotOf.erase(Fold(*keys[index]));
   delete keys[index];
   keys[index] = NULL;
   count--;
   }

StringIntHash & StringIntHash::operator = (const StringIntHash & rhs)
   {
   if (this != &rhs)
      {
      Clear();
      for (int i = 0; i < rhs.Capacity(); i++)
         if (rhs.SlotInUse(i))
            Add(*rhs.keys[i], rhs.values[i]);
      }
   return *this;
   }

//////////////////////////////////////////////////////////////////////
// StringMap and StringIntMap: sorted arrays, searched by bisection
//////////////////////////////////////////////////////////////////////

// The position of the first key that is not before s; found tells whether it equals s
static int LocateKey(const std::vector<String *> & keys, const String & s, bool & found)
   {
   int low = 0, high = (int) keys.size();

   while (low < high)
      {
      int middle = (low + high) / 2;

      if (keys[middle]->SlowCompare(s) < 0)
         low = middle + 1;
      else
         high = middle;
      }

   found = low < (int) keys.size() && keys[low]->SlowCompare(s) == 0;
   return low;
   }

// -1: no key starts with the stem, -2: more than one does
static int StemMatch(const std::vector<String *> & keys, const String & stem)
   {
   int match = -1;

   for (size_t i = 0; i < keys.size(); i++)
      if (keys[i]->SlowCompareToStem(stem) == 0)
         {
         if (match != -1) return -2;
         match = (int) i;
         }

   return match;
   }

StringMap::StringMap(int)
   {
   }

StringMap::~StringMap()
   {
   Clear();
   }

void StringMap::Clear()
   {
   for (size_t i = 0; i < keys.size(); i++)
      delete keys[i];

   keys.clear();
   objects.clear();
   }

int StringMap::Locate(const String & s, bool & found) const
   {
   return LocateKey(keys, s, found);
   }

int StringMap::Add(const String & s, void * object)
   {
   bool found;
   int  at = Locate(s, found);

   if (found)
      objects[at] = object;
   else
      {
      keys.insert(keys.begin() + at, new String(s));
      objects.insert(objects.begin() + at, object);
      }

   return at;
   }

int StringMap::Find(const String & s) const
   {
   bool found;
   int  at = Locate(s, found);

   return found ? at : -1;
   }

int StringMap::FindStem(const String & stem) const
   {
   return StemMatch(keys, stem);
   }

void StringMap::Delete(int index)
   {
   delete keys[index];
   keys.erase(keys.begin() + index);
   objects.erase(objects.begin() + index);
   }

StringIntMap::StringIntMap(int)
   {
   }

StringIntMap::~StringIntMap()
   {
   Clear();
   }

void StringIntMap::Clear()
   {
   for (size_t i = 0; i < keys.size(); i++)
      delete keys[i];

   keys.clear();
   values.clear();
   }

int StringIntMap::Locate(const String & s, bool & found) const
   {
   return LocateKey(keys, s, found);
   }

int StringIntMap::Add(const String & s, int integer)
   {
   bool found;
   int  at = Locate(s, found);

   if (found)
      values[at] = integer;
   else
      {
      keys.insert(keys.begin() + at, new String(s));
      values.insert(values.begin() + at, integer);
      }

   return at;
   }

int StringIntMap::Find(const String & s) const
   {
   bool found;
   int  at = Locate(s, found);

   return found ? at : -1;
   }

int StringIntMap::Find(const String & s, int defaultValue)
   {
   int at = Find(s);

   return at >= 0 ? at : Add(s, defaultValue);
   }

int StringIntMap::FindStem(const String & stem) const
   {
   return StemMatch(keys, stem);
   }

int StringIntMap::Integer(const String & key) const
   {
   int at = Find(key);

   return at >= 0 ? values[at] : -1;
   }

int StringIntMap::GetCount(const String & key) const
   {
   int at = Find(key);

   return at >= 0 ? values[at] : 0;
   }

int StringIntMap::IncrementCount(const String & key)
   {
   int at = Find(key);

   if (at < 0)
      {
      Add(key, 1);
      return 1;
      }

   return ++values[at];
   }

int StringIntMap::DecrementCount(const String & key)
   {
   int at = Find(key);

   if (at < 0)
      {
      Add(key, -1);
      return -1;
      }

   return --values[at];
   }

void StringIntMap::Delete(int index)
   {
   delete keys[index];
   keys.erase(keys.begin() + index);
   values.erase(values.begin() + index);
   }
