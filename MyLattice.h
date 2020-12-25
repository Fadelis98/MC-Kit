#ifndef MyLattice_H_
#define MyLattice_H_

#include <random>
#include <vector>
#include <initializer_list>

namespace lat
{
    
#pragma region 2维晶格
    static const int __2D_default_s = 100;
    static std::random_device __2D__rd__stat;
    static std::default_random_engine __2D__re__stat(__2D__rd__stat());
    static std::uniform_real_distribution<double> __2D__ur__stat;
    static std::uniform_int_distribution<int> __2D__ui__stat;
#pragma region 定义晶格类
    template <typename Type>
    class Lattice2D
    {
    private:
        int _scale;                //晶格规模
        std::vector<std::vector<Type>> _lat; //晶格数据
    public:
        Lattice2D(int s = __2D_default_s);                         //指定规模s初始化为全零向量
        Lattice2D(int s,Type val);                                //用指定的val填充晶格初始化,不可以省略s
        Lattice2D(std::initializer_list<Type>, int s = __2D_default_s); //使用参数列表初始化
        Lattice2D(const Lattice2D &);                              //复制构造函数
        Type Get(int i, int j) const;                              //周期性边界条件地选取
        Type Get(int i) const;                                     //线性索引
        Type &iGet(int i, int j);                                  //获取格点数据的引用
        Type &iGet(int i);                                         //先行索引格点数据的引用
        std::vector<std::vector<Type>> &Lat();                               //返回晶格数组的引用
        std::vector<std::vector<Type>> Data() const;                         //返回完整的晶格数据
        int Scale() const;                                         //返回当前的规模
        void Set(int i, int j, Type val);                          //设置指定格点的元素
        void SetLattice(std::vector<std::vector<Type>>);                     //用二维数组直接设置晶格数据
    };
    typedef Lattice2D<int> Ising2D;
    typedef Lattice2D<double> XY2D;
#pragma endregion
#pragma region 构造函数
    template <typename Type>
    Lattice2D<Type>::Lattice2D(int s) //默认构造函数
    {
        _scale = s;
        _lat = std::vector<std::vector<Type>>(s, std::vector<Type>(s, Type()));
    }

    template <typename Type>
    Lattice2D<Type>::Lattice2D(int s,Type val) //用指定的val填充晶格
    {
        _scale = s;
        _lat = std::vector<std::vector<Type>>(s, std::vector<Type>(s, val));
    }

    template <typename Type>
    Lattice2D<Type>::Lattice2D(std::initializer_list<Type> lst, int s) //使用参数列表初始化,每个格点将等概率地从参数列表中获取初值
    {
        _scale = s; //使用默认晶格规模
        std::vector<Type> tempvec(lst);
        int templen = lst.size();
        _lat = std::vector<std::vector<Type>>(s, std::vector<Type>(s));
        for (size_t i = 0; i < s; i++)
        {
            for (size_t j = 0; j < s; j++)
            {
                int n = __2D__ui__stat(__2D__re__stat) % templen;
                _lat[i][j] = tempvec[n];
            }
        }
    }

    template <typename Type>
    Lattice2D<Type>::Lattice2D(const Lattice2D<Type> &lat2d) //复制构造函数
    {
        _scale = lat2d.Scale();
        _lat = lat2d.Data();
    }
#pragma endregion
#pragma region 获取数据
    template <typename Type>
    Type Lattice2D<Type>::Get(int i, int j) const //周期性边界条件地获取指定格点的元素
    {
        while (i < 0)
        {
            i += _scale;
        }
        while (i >= _scale)
        {
            i -= _scale;
        }
        while (j < 0)
        {
            j += _scale;
        }
        while (j >= _scale)
        {
            j -= _scale;
        }
        return _lat[i][j];
    }

    template <typename Type>
    Type &Lattice2D<Type>::iGet(int i, int j) //周期性边界条件地获取指定格点的元素的引用
    {
        while (i < 0)
        {
            i += _scale;
        }
        while (i >= _scale)
        {
            i -= _scale;
        }
        while (j < 0)
        {
            j += _scale;
        }
        while (j >= _scale)
        {
            j -= _scale;
        }
        return _lat[i][j];
    }

    template <typename Type>
    Type Lattice2D<Type>::Get(int i) const
    {
        int indexj = i % _scale;
        int indexi = i / _scale;
        return _lat[indexi][indexj];
    }

    template <typename Type>
    Type &Lattice2D<Type>::iGet(int i)
    {
        int indexj = i % _scale;
        int indexi = i / _scale;
        return _lat[indexi][indexj];
    }

    template <typename Type>
    std::vector<std::vector<Type>> Lattice2D<Type>::Data() const //返回所有晶格数据
    {
        return _lat;
    }

    template <typename Type>
    std::vector<std::vector<Type>> &Lattice2D<Type>::Lat() //返回所有晶格数据
    {
        return _lat;
    }

    template <typename Type>
    int Lattice2D<Type>::Scale() const //返回晶格规模
    {
        return _scale;
    }
#pragma endregion
#pragma region 设置数据
    template <typename Type>
    void Lattice2D<Type>::Set(int i, int j, Type val) //设置指定格点的元素
    {
        while (i < 0)
        {
            i += _scale;
        }
        while (i >= _scale)
        {
            i -= _scale;
        }
        while (j < 0)
        {
            j += _scale;
        }
        while (j >= _scale)
        {
            j -= _scale;
        }
        _lat[i][j] = val;
    }

    template <typename Type>
    void Lattice2D<Type>::SetLattice(std::vector<std::vector<Type>> lat) //用二维数组直接设置晶格数据
    {
        _lat = lat;
    }
#pragma endregion
#pragma endregion

#pragma region 3维晶格
    static const int __3D_default_s = 10;
    static std::random_device __3D__rd__stat;
    static std::default_random_engine __3D__re__stat(__3D__rd__stat());
    static std::uniform_real_distribution<double> __3D__ur__stat;
    static std::uniform_int_distribution<int> __3D__ui__stat;
#pragma region 定义晶格类
    template <typename Type>
    class Lattice3D
    {
    private:
        int _scale;                        //晶格规模
        std::vector<std::vector<std::vector<Type>>> _lat; //晶格数据
    public:
        Lattice3D(int s = __3D_default_s);                         //指定规模s初始化为全零向量
        Lattice3D(Type val, int s);                                //用指定的val填充晶格初始化
        Lattice3D(std::initializer_list<Type>, int s = __3D_default_s); //使用参数列表初始化
        Lattice3D(const Lattice3D &);                              //复制构造函数
        Type Get(int i, int j, int k) const;                       //周期性边界条件地选取
        Type Get(int i) const;                                     //线性索引
        Type &iGet(int i, int j, int k);                           //获取格点数据的引用
        Type &iGet(int i);                                         //先行索引格点数据的引用
        std::vector<std::vector<std::vector<Type>>> &Lat();                       //返回晶格数组的引用
        std::vector<std::vector<std::vector<Type>>> Data() const;                 //返回完整的晶格数据
        int Scale() const;                                         //返回当前的规模
        void Set(int i, int j, int k, Type val);                   //设置指定格点的元素
        void SetLattice(std::vector<std::vector<std::vector<Type>>>);             //用二维数组直接设置晶格数据
    };
    typedef Lattice3D<int> Ising3D;
    typedef Lattice3D<double> XY3D;
#pragma endregion
#pragma region 构造函数
    template <typename Type>
    Lattice3D<Type>::Lattice3D(int s) //默认构造函数
    {
        _scale = s;
        _lat = std::vector<std::vector<std::vector<Type>>>(s, std::vector<std::vector<Type>>(s, std::vector<Type>(s, Type())));
    }

    template <typename Type>
    Lattice3D<Type>::Lattice3D(Type val, int s) //用指定的val填充晶格
    {
        _scale = s;
        _lat = std::vector<std::vector<std::vector<Type>>>(s, std::vector<std::vector<Type>>(s, std::vector<Type>(s, val)));
    }

    template <typename Type>
    Lattice3D<Type>::Lattice3D(std::initializer_list<Type> lst, int s) //使用参数列表初始化,每个格点将等概率地从参数列表中获取初值
    {
        _scale = s; //使用默认晶格规模
        std::vector<Type> tempvec(lst);
        int templen = lst.size();
        _lat = std::vector<std::vector<std::vector<Type>>>(s, std::vector<std::vector<Type>>(s, std::vector<Type>(s)));
        for (size_t i = 0; i < s; i++)
        {
            for (size_t j = 0; j < s; j++)
            {
                for (size_t k = 0; k < s; k++)
                {
                    int n = __3D__ui__stat(__3D__re__stat) % templen;
                    _lat[i][j][k] = tempvec[n];
                }
            }
        }
    }

    template <typename Type>
    Lattice3D<Type>::Lattice3D(const Lattice3D<Type> &lat) //复制构造函数
    {
        _scale = lat.Scale();
        _lat = lat.Data();
    }

#pragma endregion
#pragma region 获取数据
    template <typename Type>
    Type Lattice3D<Type>::Get(int i, int j, int k) const //周期性边界条件地获取指定格点的元素
    {
        while (i < 0)
        {
            i += _scale;
        }
        while (i >= _scale)
        {
            i -= _scale;
        }
        while (j < 0)
        {
            j += _scale;
        }
        while (j >= _scale)
        {
            j -= _scale;
        }
        while (k < 0)
        {
            k += _scale;
        }
        while (k >= _scale)
        {
            k -= _scale;
        }
        return _lat[i][j][k];
    }

    template <typename Type>
    Type &Lattice3D<Type>::iGet(int i, int j, int k) //周期性边界条件地获取指定格点的元素的引用
    {
        while (i < 0)
        {
            i += _scale;
        }
        while (i >= _scale)
        {
            i -= _scale;
        }
        while (j < 0)
        {
            j += _scale;
        }
        while (j >= _scale)
        {
            j -= _scale;
        }
        while (k < 0)
        {
            k += _scale;
        }
        while (k >= _scale)
        {
            k -= _scale;
        }
        return _lat[i][j][k];
    }

    template <typename Type>
    Type Lattice3D<Type>::Get(int i) const
    {
        int a = i / _scale * _scale; //层数
        int b = i % _scale * _scale; //所在层的二维线性索引
        int c = b / _scale;          //行数
        int d = b % _scale;          //列数
        int indexk = a;
        int indexj = d;
        int indexi = c;
        return _lat[indexi][indexj][indexk];
    }

    template <typename Type>
    Type &Lattice3D<Type>::iGet(int i)
    {
        int a = i / _scale * _scale; //层数
        int b = i % _scale * _scale; //所在层的二维线性索引
        int c = b / _scale;          //行数
        int d = b % _scale;          //列数
        int indexk = a;
        int indexj = d;
        int indexi = c;
        return _lat[indexi][indexj][indexk];
    }

    template <typename Type>
    std::vector<std::vector<std::vector<Type>>> Lattice3D<Type>::Data() const //返回所有晶格数据
    {
        return _lat;
    }

    template <typename Type>
    std::vector<std::vector<std::vector<Type>>> &Lattice3D<Type>::Lat() //返回所有晶格数据
    {
        return _lat;
    }

    template <typename Type>
    int Lattice3D<Type>::Scale() const //返回晶格规模
    {
        return _scale;
    }
#pragma endregion
#pragma region 设置数据
    template <typename Type>
    void Lattice3D<Type>::Set(int i, int j, int k, Type val) //设置指定格点的元素
    {
        while (i < 0)
        {
            i += _scale;
        }
        while (i >= _scale)
        {
            i -= _scale;
        }
        while (j < 0)
        {
            j += _scale;
        }
        while (j >= _scale)
        {
            j -= _scale;
        }
        while (k < 0)
        {
            k += _scale;
        }
        while (k >= _scale)
        {
            k -= _scale;
        }
        _lat[i][j][k] = val;
    }

    template <typename Type>
    void Lattice3D<Type>::SetLattice(std::vector<std::vector<std::vector<Type>>> lat) //用二维数组直接设置晶格数据
    {
        _lat = lat;
    }
#pragma endregion
#pragma endregion

} // namespace lat
#endif