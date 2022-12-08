#include "lidar_backend.h"

bool FindNearest2ID(const std_msgs::Header Headers[(WINDOW_SIZE + 1)],
                    const double tl,
                    int &id_a, int &id_b)
{
    std::vector<double> header;
    for (int i = 0; i < WINDOW_SIZE + 1; i++)
    {
        header.push_back(Headers[i].stamp.toSec());
    }
    header.push_back(tl);
    std::sort(header.begin(), header.end(), less<double>());
    // std::vector<double>::iterator it;
    auto it = std::find(header.begin(), header.end(), tl);
    if (it != header.end())
    {
        // int index = std::distance(header.begin(),it);
        int index = it - header.begin();
        id_a = index - 1;
        id_b = index;
        if (id_b > WINDOW_SIZE || id_a < 0) //|| (tl - header[id_a]) > 0.5
        {
            cout << "Not KEY LIDAR!!!" << endl;
            return false;
        }
        cout << "ID FIND!!! " << id_a << " " << id_b << endl;
        cout << "TIME CHECK!!! " << tl - header[id_a] << " " << header[id_b + 1] - tl << " " << endl;

        return true;
    }
    else
    {
        return false;
    }
}

bool FindWindowsID(const std_msgs::Header Headers[(WINDOW_SIZE + 1)],
                   const double ta, const double tb, const double tc, const double td,
                   int &id_a, int &id_b, int &id_c, int &id_d)
{
    // cout << "Windows start from " << fixed << Headers[0].stamp.toSec() << " "
    //      << "end at " << Headers[WINDOW_SIZE].stamp.toSec() << endl;
    // cout << "ta/b/c/d: " << fixed << ta << " " << tb << " " << tc << " " << td << endl;

    if ((Headers[0].stamp.toSec() > ta) || (Headers[WINDOW_SIZE].stamp.toSec() < td) || ((tb - ta) > 0.5))
    {
        cout << "Not in windows!!" << endl;
        return false;
    }

    std::vector<double> header;
    for (int i = 0; i < WINDOW_SIZE + 1; i++)
    {
        header.push_back(Headers[i].stamp.toSec());
    }

    auto ita = std::find(header.begin(), header.end(), ta);
    auto itb = std::find(header.begin(), header.end(), tb);
    auto itc = std::find(header.begin(), header.end(), tc);
    auto itd = std::find(header.begin(), header.end(), td);

    if (ita != header.end())
    {
        id_a = std::distance(header.begin(), ita);
    }
    if (itb != header.end())
    {
        id_b = std::distance(header.begin(), itb);
    }
    if (itc != header.end())
    {
        id_c = std::distance(header.begin(), itc);
    }
    if (itd != header.end())
    {
        id_d = std::distance(header.begin(), itd);
    }
    if (id_b == id_c)
    {
        id_a--;
        id_b--;
    }

    if (id_b > id_a && id_d > id_c && id_a >= 0 && id_a != id_c)
    {
        return true;
    }
    else
    {
        return false;
    }
}
