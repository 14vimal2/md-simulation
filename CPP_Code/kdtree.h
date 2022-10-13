#include <iostream>

using namespace std;

/************************* Median of Median *****************************
 * adapted from @link https://en.wikipedia.org/wiki/Median_of_medians
 *
 * **********************************************************************/
template <typename _type_>
void kswap(_type_ *a, _type_ *b, uint8_t k)
{
    for (uint8_t i = 0; i < k; i++)
    {
        swap(a[i], b[i]);
    }
}

template <typename _dec_type_, typename _int_type_>
_int_type_ pivot(_dec_type_ *arr, _int_type_ left, _int_type_ right, uint8_t axis, uint8_t k);

template <typename _dec_type_, typename _int_type_>
_int_type_ partition5(_dec_type_ *arr, _int_type_ left, _int_type_ right, uint8_t axis, uint8_t k)
{
    _int_type_ i = left + k, j;
    while (i <= right)
    {
        j = i;
        while (j > left && arr[j - k + axis] > arr[j + axis])
        {
            kswap(arr + j - k, arr + j, k);
            j -= k;
        }
        i += k;
    }
    return ((left + right) / 2 / k) * k;
}

template <typename _dec_type_, typename _int_type_>
_int_type_ partition(_dec_type_ *arr, _int_type_ left, _int_type_ right, _int_type_ pivotIndex, _int_type_ n, uint8_t axis, uint8_t k)
{
    _dec_type_ pivotValue = arr[pivotIndex + axis];
    kswap(arr + pivotIndex, arr + right, k);
    _int_type_ storeIndex = left;

    for (_int_type_ i = left; i < right; i += k)
    {
        if (arr[i + axis] < pivotValue)
        {
            kswap(arr + storeIndex, arr + i, k);
            storeIndex += k;
        }
    }

    _int_type_ storeIndexEq = storeIndex;

    for (_int_type_ i = storeIndex; i < right; i += k)
    {
        if (arr[i + axis] == pivotValue)
        {
            kswap(arr + storeIndexEq, arr + i, k);
            storeIndexEq += k;
        }
    }

    kswap(arr + right, arr + storeIndexEq, k);

    if (n < storeIndex)
    {
        return storeIndex;
    }

    if (n <= storeIndexEq)
    {
        return n;
    }

    return storeIndexEq;
}

template <typename _dec_type_, typename _int_type_>
_int_type_ select(_dec_type_ *arr, _int_type_ left, _int_type_ right, _int_type_ n, uint8_t axis, uint8_t k)
{
    _int_type_ pivotIndex;
    while (left != right)
    {
        pivotIndex = pivot(arr, left, right, axis, k);
        pivotIndex = partition(arr, left, right, pivotIndex, n, axis, k);
        if (n == pivotIndex)
        {
            return n;
        }
        else if (n < pivotIndex)
        {
            right = pivotIndex - k;
        }
        else
        {
            left = pivotIndex + k;
        }
    }

    return left;
}

template <typename _dec_type_, typename _int_type_>
_int_type_ pivot(_dec_type_ *arr, _int_type_ left, _int_type_ right, uint8_t axis, uint8_t k)
{
    if (right - left < 5 * k)
    {
        return partition5(arr, left, right, axis, k);
    }
    _int_type_ subRight, median5;
    for (_int_type_ i = left; i < right; i += 5 * k)
    {
        subRight = i + 4 * k;
        if (subRight > right)
        {
            subRight = right;
        }

        median5 = partition5(arr, i, subRight, axis, k);
        kswap(arr + median5, arr + left + ((i - left) / 5 / k) * k, k);
    }

    _int_type_ mid = ((right - left) / 10 / k) * k + left + k;
    return select(arr, left, left + ((right - left) / 5 / k) * k, mid, axis, k);
};

