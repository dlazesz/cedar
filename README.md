Cedar
=====
This is a C++ implementation of efficiently-updatable double-array trie, developed by Naoki Yoshinaga at Kitsuregawa, Toyoda Lab., IIS, University of Tokyo.

If you make use of cedar for research or commercial purposes, the reference will be:

  N. Yoshinaga and M. Kitsuregawa. A Self-adaptive Classifier for Efficient Text-stream Processing. Proc. COLING 2014, 2014.

Modifications
-------------

All modifications originate form Cedar (June 24th 2014):
- Cedar had a limitation on keys in favour of memory usage. **This modified version is able to handle more than 10^9 keys, but has more memory usage.** The trie must be declared with `long` type insted of `int` for the reduced trie to work!
- Build script has changed from **automake** to **cmake**.
- Efficient reset method to reset the trie without reallocating the memory (from: https://github.com/KrishnaPG/cedar)
- Additional CommonPrefixSearch() based on sentinel (without the need for computing the string length) (from: https://github.com/KrishnaPG/cedar)
- Add option to set memory upperbound and change behaviour to try to allocate less memory if the allocation of double amount is failed.

**Keys with `\00` in them and zero length keys still not supported!**

License
======
GNU GPLv2, LGPLv2.1, and BSD (see License file)

Additionally the original code is also vailable under GPL/LGPL terms at the website location given below.

Author
======
Naoki Yoshinaga at Kitsuregawa, Toyoda Lab., IIS, University of Tokyo

Website: http://www.tkl.iis.u-tokyo.ac.jp/~ynaga/cedar

Setup
======

1. $ mkdir build && cd build
2. $ cmake .. && make && sudo make install


For using as a library the only needed files are:

- cedar.h
- cedarpp.h

There are standalone tools:

- mkcedar: Create tree from text file
- cedar: Interactive demo on different search functions
- simple.cc is a simple demo on usage (not installed with make install)

For detailed API reference visit the website: http://www.tkl.iis.u-tokyo.ac.jp/~ynaga/cedar
