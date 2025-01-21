echo "";
echo "";

echo "~~ clearing old results... ~~";

echo "";
echo "---->  rm ./renders/*.png";
echo "---->  rm ./res/*.dat";
rm ./renders/*.png;
rm ./res/*.dat;
echo "";

echo "~~ compilling... ~~";

echo "";
echo "---->  g++ main.cpp -o main";
echo "";
g++ main.cpp -O2 -o main;

echo "~~ compilled !, now running  ~~";

echo "";
echo "----> ./main";
echo "";
./main;

echo "~~ values computed, now plotting ~~";

echo "";
echo "----> python3 modules/plotting.py";
echo "";
python3 modules/plotting.py
