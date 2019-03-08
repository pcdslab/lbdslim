/*
  The MIT License (MIT)

Copyright (c) 2016 Eduard López

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef KEYVAL_H_
#define KEYVAL_H_

#include <omp.h>
#include "common.h"

/*
 * FUNCTION: KeyVal_Serial_Internal
 *
 * DESCRIPTION: Internal Serial Key Value Sort
 *
 * INPUT:
 * @key  : Key array
 * @val  : Value array
 * @left : Start index of arrays
 * @right: End index of arrays
 *
 * OUTPUT: none
 */
template<class K, class V>
VOID KeyVal_Serial_Internal(K *key, V *val, INT low, INT high)
{
    INT i = low;
    INT j = high;
    K pivotkey = key[(i + j) / 2];
    K tmpkey;
    V tmpval;

    while (i <= j)
    {
        while (key[i] < pivotkey)
            i++;
        while (key[j] > pivotkey)
            j--;

        if (i <= j)
        {
            tmpkey = key[i];
            tmpval = val[i];

            key[i] = key[j];
            val[i] = val[j];

            key[j] = tmpkey;
            val[j] = tmpval;

            i++;
            j--;
        }
    }

    if (j > low)
        KeyVal_Serial_Internal(key, val, low, j);
    if (i < high)
        KeyVal_Serial_Internal(key, val, i, high);
}

/*
 * FUNCTION: KeyVal_Serial
 *
 * DESCRIPTION: Internal Parallel Key Value Sort
 *
 * INPUT:
 * @key     : Key array
 * @val     : Value array
 * @lenArray: Length of arrays
 *
 * OUTPUT: none
 */
template<class K, class V>
VOID KeyVal_Serial(K *key, V *val, UINT lenArray)
{
    KeyVal_Serial_Internal<K, V>(key, val, 0, lenArray-1);
}


#endif /* KEYVAL_H_ */
