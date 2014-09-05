開度作成プログラム
=====
地上開度、地下開度、尾根谷度を作成するためのプログラムです。
gdalとOpenMPを利用しています。

#使い方
```
kaido_mp.exe input.tif output.tif radius type_no thread_count

radius 集計半径
type_no 地上開度:1 地下開度:2 尾根谷度:3
thread_count マルチスレッド数
```


例 半径250mの地上開度を4スレッドで計算

```
kaido_mp.exe demUTM.tif tijyoukaido.tif 250 1 4
```

#参考
- https://www.jstage.jst.go.jp/article/jsprs1975/38/4/38_4_26/_article/-char/ja/
- http://archive.sokugikyo.or.jp/pdf/apa96_2008_03/APA963.pdf

#ライセンス&注意

- このプログラムはMIT Licenseです。
- プログラムにバグあっても責任を負いません。
