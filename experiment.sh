#!/bin/bash

# 定义要执行的测试文件列表（按顺序排列）
test_files=(
    "test2groups"
    "test3groups"
    "test4groups"
    "test5groups"
    "test6groups"
    "test7groups"
    "test8groups"
)

# 遍历所有测试文件，依次执行
for file in "${test_files[@]}"; do
    # 检查文件是否存在且可执行
    if [ -x "./build/bin/examples/binfhe/$file" ]; then
        echo -e "\n====================================="
        echo "开始执行：$file"
        echo "====================================="
        
        # 执行文件（$file 前加 ./ 确保执行当前目录下的文件）
        ./build/bin/examples/binfhe/"$file"
        
        # 捕获执行结果，提示成功/失败
        if [ $? -eq 0 ]; then
            echo -e "\n✅ $file 执行成功！"
        else
            echo -e "\n❌ $file 执行失败！退出码：$?"
        fi
    else
        # 文件不存在或不可执行时的提示
        echo -e "\n⚠️  跳过：$file 不存在或无执行权限"
    fi
done

echo -e "\n====================================="
echo "所有测试文件执行完毕！"
echo "====================================="