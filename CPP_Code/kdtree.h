
#include <iostream>
#ifndef _INCL_KDTREE_
#define _INCL_KDTREE_
using namespace std;

/************************* Median of Medians *****************************
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

//************************* End of Median of Medians **************************//

template <typename _dec_type_, typename _int_type_>
class kdtree
{
private:
    class Node
    {
    private:
        _int_type_ ind;
        Node *left, *right;

    public:
        Node(_int_type_ x)
        {
            ind = x;
            left = nullptr;
            right = nullptr;
        };

        bool isSame(_dec_type_ *point, _dec_type_ *data_, _int_type_ k)
        {
            for (_int_type_ i = 0; i < k; i++)
            {
                if (point[i] != data_[ind + i])
                {
                    return false;
                }
            }
            return true;
        }

        _dec_type_ distance2(_dec_type_ *point, _dec_type_ *data_, _int_type_ k)
        {
            _dec_type_ dist = 0;
            for (_int_type_ i = 0; i < k; i++)
            {
                dist += (data_[ind + i] - point[i]) * (data_[ind + i] - point[i]);
            }
            return dist;
        }

        // returns true the given point is in the left of the node
        bool inLeft(_dec_type_ *point, _dec_type_ *data_, _int_type_ axis)
        {
            return point[axis] < data_[ind + axis];
        }

        bool inRange(_dec_type_ *point, _dec_type_ *data_, _int_type_ axis, _dec_type_ radius)
        {
            return (point[axis] - radius < data_[ind + axis]) && (point[axis] + radius > data_[ind + axis]);
        }

        void searchNeighbours(_dec_type_ *point, _dec_type_ *forceVal, _dec_type_ *data_, _dec_type_ radius, _int_type_ k, _int_type_ depth = 0)
        {
            // cout << "points searched " << endl;
            // for (int i = 0; i < k; i++)
            // {
            //     cout << point[i] << " ";
            // }
            // cout << endl;

            _int_type_ axis = depth % k;
            if (!isSame(point, data_, k))
            {
                if (inRange(point, data_, axis, radius))
                {
                    _dec_type_ r2 = distance2(point, data_, k);
                    if (r2 < radius * radius)
                    {
                        // for (uint8_t i = 0; i < k; i++)
                        // {
                        //     cout << data_[ind + i] << " ";
                        // }
                        // cout << endl;

                        _dec_type_ r2i = 1.0 / r2, r6i;
                        r6i = r2i * r2i * r2i;
                        _dec_type_ ff = 48.0 * r2i * r6i * (r6i - 0.5);

                        for (uint8_t i = 0; i < k; i++)
                        {
                            forceVal[i] += ff * (point[i] - data_[ind + i]);
                        }
                    }
                    if (left)
                    {
                        left->searchNeighbours(point, forceVal, data_, radius, k, ++depth);
                    }
                    if (right)
                    {
                        right->searchNeighbours(point, forceVal, data_, radius, k, ++depth);
                    }
                    return;
                }
                else if (inLeft(point, data_, axis))
                {
                    if (left)
                    {
                        left->searchNeighbours(point, forceVal, data_, radius, k, ++depth);
                    }
                }
                else if (right)
                {
                    right->searchNeighbours(point, forceVal, data_, radius, k, ++depth);
                }
            }
        }

        void setleft(Node *L) { left = L; }
        void setright(Node *R) { right = R; }
        void deleteNode()
        {
            if (left)
            {
                left->deleteNode();
            }
            delete left;
            if (right)
            {
                right->deleteNode();
            }
            delete right;
        }

        void displayNode(_dec_type_ *arr, _int_type_ k)
        {
            cout << "((";
            cout << arr[ind];
            for (uint8_t i = 1; i < k; i++)
            {
                cout << ", " << arr[ind + i];
            }

            cout << "), ";
            if (left)
            {
                left->displayNode(arr, k);
            }
            cout << ", ";
            if (right)
            {
                right->displayNode(arr, k);
            }
            cout << " )";
        }
    };
    _dec_type_ *data;
    _int_type_ k;
    _int_type_ N;
    Node *root;
    _dec_type_ boxsize;

public:
    kdtree(_dec_type_ *arr, _int_type_ N_, _int_type_ k_, _dec_type_ boxsize_)
    {
        if (N_ * k_)
        {
            data = new _dec_type_[N_ * k_];
            k = k_;
            N = N_;
            boxsize = boxsize_;
            copy(arr, arr + N * k, data);
            root = buildtree(0, (N - 1) * k, 0);
        }
        else
        {
            root = nullptr;
        }
    };

    void refresh(_dec_type_ *arr)
    {

        if (root)
        {
            root->deleteNode();
        }
        delete root;

        copy(arr, arr + N * k, data);
        root = buildtree(0, (N - 1) * k, 0);
    }

    Node *buildtree(_int_type_ l, _int_type_ r, _int_type_ depth)
    {
        if (l > r)
        {
            return nullptr;
        }

        if (l == r)
        {
            return new Node(l);
        }
        _int_type_ axis = depth % k;
        _int_type_ medianIndex = select(data, l, r, ((l + r + k) / 2 / k) * k, axis, k);
        // cout << "median " << medianIndex << endl;


        // displayTree();

        while (medianIndex > k &&  medianIndex > l &&  data[medianIndex + axis] == data[medianIndex + axis - k])
        {
            medianIndex -= k;
        }
        Node *temp = new Node(medianIndex);
        temp->setleft(buildtree(l, medianIndex - k, ++depth));
        temp->setright(buildtree(medianIndex + k, r, ++depth));
        return temp;
    }

    void findNeighboursInRange(_dec_type_ *point, _dec_type_ *forceval, _dec_type_ radius)
    {
        if (root)
        {
            root->searchNeighbours(point, forceval, data, radius, k);
            for (uint8_t i = 0; i < k; i++)
            {
                if (point[i] < radius || boxsize - point[i] < radius)
                {
                    _dec_type_ temp = point[i];
                    point[i] = point[i] < radius ? point[i] + boxsize : point[i] - boxsize;
                    root->searchNeighbours(point, forceval, data, radius, k);
                    point[i] = temp;
                }
            }
        }
    }
    void displayTree()
    {
        cout << endl;
        if (root)
        {
            root->displayNode(data, k);
        }
    }
    void displayData()
    {
        cout << "N : " << N << endl;
        cout << "k : " << k << endl;
        for (int i = 0; i < N * k; i++)
        {
            cout << data[i] << " ";
        }
    }
    ~kdtree()
    {
        if (root)
        {
            root->deleteNode();
        }
        delete root;
        delete[] data;
    }
};

#endif