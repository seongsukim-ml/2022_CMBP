        Standard result save
        string result = to_string(cBin) + "," + to_string(model.TV[cBin]) + "," + to_string(model.MV[cBin]) + "," + to_string(model.CV[cBin]) + ",";
        result = result + to_string(model.res[0]) + "," + to_string(model.res[1]) + "," + to_string(model.res[2]) + ",";
        result = result + to_string(model.res[3]) + "," + to_string(model.res[4]) + "," + to_string(model.BV[cBin]) + ",";
        result = result + to_string(MMerr) + "," + to_string(CCerr) + "," + to_string(BBerr) + ",";
        result = result + to_string(MM2err) + "," + to_string(MM4err) + "\n";