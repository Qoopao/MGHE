#!/bin/bash

total_runs=10
#expected="Result of encrypted computation of ( 1 NAND 1 ) = 0"
expected="---------第一个NAND门计算完成-------------
Result of encrypted computation of ( ct2:0 NAND ct3:0 ) = 1
---------第二个NAND门计算完成-------------
Result of encrypted computation of ct1:1 NAND ( ct2:0 NAND ct3:0 ) = 0
---------第三个NAND门计算完成-------------
Result of encrypted computation of ct1:1 NAND ( ct2:0 NAND ct3:0 ) NAND ct4:1 = 1"

#删除user目录下的rrr.log文件
rm -f /vscode/myProgram/RRRREALL/rrr.log

echo "运行 $total_runs 次 ./MG-boolean-xzddf 测试..."

for ((i=1; i<=total_runs; i++)); do
    result=$(./build/bin/examples/binfhe/MG-boolean-xzddf 2>&1 | tail -n 6)
    
    if [ "$result" != "$expected" ]; then
        echo "❌ 第 $i 次失败: $result"
        echo "程序停止"
        exit 1
    fi
    
    # 显示进度
    echo "已完成 $i/$total_runs 次测试"
done

echo "🎉 所有 $total_runs 次测试均成功！"