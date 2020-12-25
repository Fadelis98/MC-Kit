#include "anyMc.h"
#include <fstream>
using namespace std;
using namespace lat;

namespace lat
{
    static std::random_device __xymodel_rd;
    static std::default_random_engine __xychoose_e(__xymodel_rd());
    template <>
    double &chooseRvec(double &numOfpi)
    {
        static uniform_real_distribution<double> __xychoose_ur(0, 2);
        numOfpi = __xychoose_ur(__xychoose_e);
        return numOfpi;
    }
    template <>
    double reflect(const double &orign, const double &rvec)
    {
        double val = rvec - (orign - rvec);
        if (val >= 2)
        {
            return val - 2;
        }
        else if (val < 0)
        {
            return val + 2;
        }
        return val;
    }
    template <>
    double caculatePairEnergy(const double &sigma1, const double &sigma2)
    {
        return -cos((sigma1 - sigma2) * M_PI);
    }
    template <>
    double caculateMagnitude(const Lattice2D<double> &model)
    {
        int scale = model.Scale();
        double MagX, MagY = 0;
        for (size_t x = 0; x < scale; x++)
        {
            for (size_t y = 0; y < scale; y++)
            {
                double deg = model.Get(x, y) * M_PI;
                MagX += cos(deg);
                MagY += sin(deg);
            }
        }
        return sqrt(MagX * MagX + MagY * MagY);
    }
} // namespace lat

vector<double> TList = linspace(0, 3, 6); //需要演算的温度列表
vector<int> ScaleList = {20, 50};         //指定对哪些规模的晶格进行模拟和采样
int main()
{
    size_t numOfdataPoints = 2000; //采样点的数量
    size_t preheat_Wolff = 0;   //Wolff算法预热模型需要的次数
    size_t preheat_MP = 0;     //Metropolis算法预热模型需要的次数
    size_t intervalWolff = 10;     //wolff算法的采样间隔
    size_t intervalMP = 100;       //Metropolis算法的采样间隔

    for (auto scale : ScaleList)
    {
        for (auto temperature : TList)
        {
            double beta = 1 / temperature;
            //MetroPolis算法模拟的段落
            Lattice2D<int> MPmodel(scale);
            ofstream dataofMP("MPmag" + to_string(temperature) + ".csv");
            for (size_t i = 0; i < preheat_MP; i++) //预热模型
            {
                MetroPolisProcess(MPmodel, beta);
            }
            for (size_t i = 0; i < numOfdataPoints; i++) //采样
            {
                for (size_t j = 0; j < intervalMP; j++) //采样间隔
                {
                    MetroPolisProcess(MPmodel, beta);
                }
                double mag = caculateMagnitude(MPmodel);
                dataofMP << mag << endl;
            }
            dataofMP.close();
            //Wolff算法模拟的段落
            Lattice2D<int> Wolffmodel(scale);
            ofstream dataofWolff("Wolffmag" + to_string(temperature) + ".csv");
            for (size_t i = 0; i < preheat_Wolff; i++) //预热模型
            {
                MetroPolisProcess(Wolffmodel, beta);
            }
            for (size_t i = 0; i < numOfdataPoints; i++) //采样
            {
                for (size_t j = 0; j < intervalWolff; j++) //采样间隔
                {
                    MetroPolisProcess(Wolffmodel, beta);
                }
                double mag = caculateMagnitude(Wolffmodel);
                dataofWolff << mag << endl;
            }
            dataofMP.close();
        }
    }
}