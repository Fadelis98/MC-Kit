#ifndef ANYMC_H_
#define ANYMC_H_

#include "MyLattice.h"
#include <cmath>

namespace lat
{
    static std::random_device __AnyMc_rd;
    static std::default_random_engine __MP_e(__AnyMc_rd());
    static std::default_random_engine __wolff_e(__AnyMc_rd());
    static std::default_random_engine __wolff_ei(__AnyMc_rd());

#pragma region 声明坐标结构
    struct position2d
    {
        int x;
        int y;
    };
    struct position3d
    {
        int x;
        int y;
        int z;
    };
#pragma endregion
#pragma region 声明函数

    #pragma region 辅助函数
        //生成等间距列表
        std::vector<double> linspace(double a, double b, int num);
        //数组求模取平均
        double absmean(const std::vector<double> &);
        //数组求平均
        double mean(const std::vector<double> &vec);

    #pragma endregion

    #pragma region metropolis算法相关函数
        template <typename Type>
        bool MPAcceptOrNot(const position2d &p, const Lattice2D<Type> &model, const double beta, const Type &rvec);
        template <typename Type>
        void MetroPolisProcess(Lattice2D<Type> &model, const double beta);
    #pragma endregion

    #pragma region wolff算法相关函数
        //计算是否接受
        template <typename Type>
        bool WolffAcceptOrNot(const position2d &p, const position2d &lastp, const Lattice2D<Type> &model, const double beta, const Type &rvec);
        //递归地扩充集团
        template <typename Type>
        void recursiveExpand(const position2d &p, const position2d &lastp, Lattice2D<Type> &model, Lattice2D<int> &label, int &count, const double beta, const Type &rvec);

        template <typename Type>
        int WolffProcess(Lattice2D<Type> &model, const double beta);
    #pragma endregion

    #pragma region 具体系统相关的函数，需要用户自定义
        //定义反射矢量选取方法
        template <typename Type>
        Type& chooseRvec(Type &); //用户自定义
        //定义反射方法
        template <typename Type>
        Type reflect(const Type &orign, const Type &rvec); //用户自定义
        //定义自旋对能量
        template <typename Type>
        double caculatePairEnergy(const Type &sigma1, const Type &sigma2); //用户自定义
        //定义总磁矩
        template <typename Type>
        Type caculateMagnitude(const Lattice2D<Type> &model);
    #pragma endregion

    #pragma region 计算物理量
        //计算总能量
        template <typename Type>
        double caculateEnergy(const Lattice2D<Type> &model);
        //计算自旋能量
        template <typename Type>
        double caculateSpinEnergy(position2d p, const Lattice2D<Type> &model);
        //给定反射参数计算试探自旋能量
        template <typename Type>
        double caculateSpinEnergy(position2d p, const Lattice2D<Type> &model, Type rvec);

        //计算binder ratio
        double caculateBiner(const std::vector<double> &Magnitude, int begin = 0, int end = -1);
    #pragma endregion

#pragma endregion 声明函数

#pragma region 定义函数

    #pragma region 预定义具体模型相关函数
        template <typename Type>
        Type& chooseRvec(Type &obj)
        {
            obj = Type();
            return obj;
        }

        template <typename Type>
        Type reflect(const Type &orign, const Type &rvec)
        {
            return orign - 2 * (orign * rvec) * rvec;
        }

        template <typename Type>
        double caculatePairEnergy(const Type &sigma1, const Type &sigma2)
        {
            return -sigma1 * sigma2;
        }

        template <typename Type>
        Type caculateMagnitude(const Lattice2D<Type> &model)
        {
            int scale = model.Scale();
            Type Mag;
            for (size_t x = 0; x < scale; x++)
            {
                for (size_t y = 0; y < scale; y++)
                {
                    Mag += model.Get(x, y);
                }
            }
            return Mag;
        }

    #pragma endregion

    #pragma region 辅助函数
        std::vector<double> linspace(double a, double b, int num)
        {
            double space = (b - a) / (num - 1);
            std::vector<double> vec(num);
            for (size_t i = 0; i < num; i++)
            {
                vec[i] = (a + i * space);
            }
            return vec;
        }

        double absmean(const std::vector<double> &vec)
        {
            double mean = 0;
            for (auto i : vec)
            {
                mean += abs(i);
            }
            mean /= vec.size();
            return mean;
        }

        double mean(const std::vector<double> &vec)
        {
            double mean = 0;
            for (auto i : vec)
            {
                mean += i;
            }
            mean /= vec.size();
            return mean;
        }

        double caculateBiner(const std::vector<double> &MeanMagnitude, int begin, int end)
        {
            if (end == -1)
            {
                end = MeanMagnitude.size();
            }
            double m1 = 0;
            long double m2 = 0;
            long double m4 = 0;
            for (auto i = MeanMagnitude.begin() + begin; i < MeanMagnitude.begin() + end; i++)
            {
                double m = *i;
                double mm = m * m;
                m1 += m;
                m2 += mm;
                m4 += mm * mm;
            }
            int num = end - begin;
            m1 /= num;
            m2 /= num;
            m4 /= num;
            return 1.5 - m4 / (m2 * m2 * 2);
        }
    #pragma endregion

    #pragma region Wolff算法相关

        template <typename Type>
        bool WolffAcceptOrNot(const position2d &p, const position2d &lastp, const Lattice2D<Type> &model, const double beta, const Type &rvec)
        {
            if (std::isinf(beta))
            {
                return true;
            }
            Type newone = model.Get(p.x, p.y);
            Type lastone = model.Get(lastp.x, lastp.y);
            double OldEnergy = caculatePairEnergy(newone, lastone);
            double NewEnergy = caculatePairEnergy(reflect(newone, rvec), lastone);
            double delta = -NewEnergy + OldEnergy;
            double probability = exp(-beta * delta);
            static std::uniform_real_distribution<double> __wolff_ur;
            double randp = __wolff_ur(__wolff_e);
            return (randp > probability);
        }

        template <typename Type>
        void recursiveExpand(const position2d &p, const position2d &lastp, Lattice2D<Type> &model, Lattice2D<bool> &label, int &count, const double beta, const Type &rvec)
        {
            bool stat = label.Get(p.x, p.y);
            if (!stat) //对本次扩张中未访问过的节点
            {
                if (count == 0 || WolffAcceptOrNot(p, lastp, model, beta, rvec))
                {
                    //接受节点加入
                    label.Set(p.x, p.y, true);
                    count++;
                    //翻转节点
                    model.Set(p.x, p.y, reflect(model.Get(p.x, p.y), rvec));
                    //递归地向最近邻节点扩张
                    recursiveExpand(position2d{p.x + 1, p.y}, position2d{p.x, p.y}, model, label, count, beta, rvec);
                    recursiveExpand(position2d{p.x - 1, p.y}, position2d{p.x, p.y}, model, label, count, beta, rvec);
                    recursiveExpand(position2d{p.x, p.y + 1}, position2d{p.x, p.y}, model, label, count, beta, rvec);
                    recursiveExpand(position2d{p.x, p.y - 1}, position2d{p.x, p.y}, model, label, count, beta, rvec);
                }
                else //拒绝加入直接返回
                {
                    return;
                }
            }
            else //遇到已访问的节点则返回
            {
                return;
            }
        }
        template <typename Type>
        int WolffProcess(Lattice2D<Type> &model, const double beta)
        {
            int scale = model.Scale();
            //选择随机反射矢量
            Type rvec;
            chooseRvec(rvec);
            //选择随机种子

            std::uniform_int_distribution<int> __wolff_uui(0, scale - 1);
            position2d seed{__wolff_uui(__wolff_ei), __wolff_uui(__wolff_ei)};
            //初始化晶格标记
            Lattice2D<bool> label(scale, false);
            //翻转选取结点
            int count = 0;
            recursiveExpand(seed, seed, model, label, count, beta, rvec);
            return count;
        }

    #pragma endregion

    #pragma region MetroPolis算法相关
        template <typename Type>
        bool MPAcceptOrNot(const position2d &p, const Lattice2D<Type> &model, const double beta, const Type &rvec)
        {
            if (std::isinf(beta))
            {
                return false;
            }
            double OldEnergy = caculateSpinEnergy(p, model);
            double NewEnergy = caculateSpinEnergy(p, model, rvec);
            double delta = NewEnergy - OldEnergy;
            double probability = exp(-beta * delta);
            static std::uniform_real_distribution<double> __MP_ur;
            double randp = __MP_ur(__MP_e);
            return (randp < probability);
        }
        template <typename Type>
        void MetroPolisProcess(Lattice2D<Type> &model, const double beta)
        {
            int scale = model.Scale();
            std::uniform_int_distribution<int> __MP_ui(0, scale - 1);
            position2d seed{__MP_ui(__MP_e), __MP_ui(__MP_e)};
            Type rvec;
            chooseRvec(rvec);
            if (MPAcceptOrNot(seed, model, beta, rvec))
            {
                model.Set(seed.x, seed.y, reflect(model.Get(seed.x, seed.y), rvec));
            }
        }
    #pragma endregion

    #pragma region 计算物理量
        template <typename Type>
        double caculateEnergy(const Lattice2D<Type> &model)
        {
            int scale = model.Scale();
            double Energy = 0;
            for (size_t x = 0; x < scale; x++)
            {
                for (size_t y = 0; y < scale; y++)
                {
                    Energy += caculateSpinEnergy(position2d{x, y}, model);
                }
            }
            return Energy / 2;
        }

        template <typename Type>
        double caculateSpinEnergy(position2d p, const Lattice2D<Type> &model)
        {
            Type orign = model.Get(p.x, p.y);
            double energy = (caculatePairEnergy(orign, model.Get(p.x - 1, p.y)) +
                            caculatePairEnergy(orign, model.Get(p.x + 1, p.y)) +
                            caculatePairEnergy(orign, model.Get(p.x, p.y + 1)) +
                            caculatePairEnergy(orign, model.Get(p.x, p.y - 1)));
            return energy;
        }

        template <typename Type>
        double caculateSpinEnergy(position2d p, const Lattice2D<Type> &model, Type rvec)
        {
            Type changed = reflect(model.Get(p.x, p.y), rvec);
            double energy = (caculatePairEnergy(changed, model.Get(p.x - 1, p.y)) +
                            caculatePairEnergy(changed, model.Get(p.x + 1, p.y)) +
                            caculatePairEnergy(changed, model.Get(p.x, p.y + 1)) +
                            caculatePairEnergy(changed, model.Get(p.x, p.y - 1)));
            return energy;
        }

    #pragma endregion

#pragma endregion

} // namespace lat
#endif