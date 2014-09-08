#コンパイル方法

##windowsの場合
1. VC2010をインストール　http://go.microsoft.com/fwlink/?LinkId=190491
2. OpenMP環境インストール() http://www.microsoft.com/en-us/download/confirmation.aspx?id=11310
3. gdalをインストール
4. OpenMPのためのパスを設定
```
Set INCLUDE=C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\INCLUDE;%INCLUDE%
Set LIB=C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\LIB;%LIB%
```
5. Visual Studio コマンドプロンプト(2010)でコンパイル
```
cl /EHsc /openmp -IC:\warmerda\bld\include C:\warmerda\bld\lib\gdal_i.lib kaido_mp.cpp
```

##macの場合

1. homebrewをインストール
2. brew install gdal
3. コンパイル
```
g++ -I/usr/local/Cellar/gdal/1.10.1/include/ -L/usr/local/Cellar/gdal/1.10.1/lib kaido.cpp /usr/local/Cellar/gdal/1.10.1/lib/libgdal.dylib -fopenmp -o kaido
```