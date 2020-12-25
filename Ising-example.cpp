#include "anyMc.h"
#include <fstream>
using namespace std;
using namespace lat;

namespace lat
{
    template <>
    int& chooseRvec(int &numOfpi)
    {
        return numOfpi;
    }
    template <>
    int reflect(const int &orign, const int &rvec)
    {
        return orign*-1;
    }
    template <>
    double caculatePairEnergy(const int &sigma1, const int &sigma2)
    {
        return -sigma1*sigma2;
    }
} // namespace lat

vector<double> TList = linspace(0,3,6);//需要演算的温度列表
vector<int> ScaleList = {20,50};//指定对哪些规模的晶格进行模拟和采样
int main()
{
    int numOfdataPoints = 2000;//采样点的数量
    int preheat_iter = 0;      //预热模型需要的次数
    int intervalWolff = 10;    //wolff算法的采样间隔
    int intervalMP = 100;      //Metropolis算法的采样间隔
    for (auto scale : ScaleList)
    {
        for (auto temperature : TList)
        {
            double beta = 1 / temperature;
            //MetroPolis算法模拟的段落
            XY2D MPmodel(scale);
            ofstream dataofMP("MPmag" + to_string(temperature) + ".csv");
            for (size_t i = 0; i < preheat_iter; i++) //预热模型
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
            XY2D Wolffmodel(scale);
            ofstream dataofWolff("Wolffmag" + to_string(temperature) + ".csv");
            for (size_t i = 0; i < preheat_iter; i++) //预热模型
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