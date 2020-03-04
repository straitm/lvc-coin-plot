#!/bin/bash

mfdndligo()
{
  e=$1
  mass=$2
  ligopass2/ligopass2.sh ligoresults-$e/20*fardet-ligo.hadded.root | \
    grep snprob$mass | awk '{print $2, $3}' > a$e &
  ligopass2/ligopass2.sh ligoresults-$e/20*neardet-ligo.hadded.root | \
    grep snprob$mass | awk '{print $3}' > b$e &

  wait

  paste a$e b$e > comb$e

  printf '%4s fdndligo   %11s %7s of 10kpc, %7s kpc %s %s\n' \
    $mass \
    $e \
    $(ligopass2/find90 $e $mass < comb$e | awk '{print $4, $9, $11, $12}')
}

fdndligo()
{
  mfdndligo $1 27
  mfdndligo $1 9.6
}

mt02ndligo()
{
  e=$1
  mass=$2
  ligopass2/ligopass2.sh ligoresults-$e/20*fardet-t02.hadded.root | \
    grep snprob$mass | awk '{print $2, $3}' > a$e &
  ligopass2/ligopass2.sh ligoresults-$e/20*neardet-ligo.hadded.root | \
    grep snprob$mass | awk '{print $3}' > b$e &

  wait

  paste a$e b$e > comb$e

  printf '%4s t02+NDligo %11s %7s of 10kpc, %7s kpc %s %s\n' \
    $mass \
    $e \
    $(ligopass2/find90 $e $mass < comb$e | awk '{print $4, $9, $11, $12}')
}

t02ndligo()
{
  mt02ndligo $1 27
  mt02ndligo $1 9.6
}

mjustt02()
{
  e=$1
  mass=$2
  printf '%4s FD-t02     %11s %7s of 10kpc, %7s kpc\n' \
    $mass \
    $e \
    $(ligopass2/ligopass2.sh ligoresults-$e/20*fardet-t02.hadded.root | \
      grep "${mass}SN-like flux" | \
      awk '{print $4, $9}')
}

justt02()
{
  mjustt02 $1 27
  mjustt02 $1 9.6
}

ligopass2/ligopass2.sh # no args means compile

t02ndligo S190426c
t02ndligo S190513bm
fdndligo  S190521g
fdndligo  S190602aq
fdndligo  S190630ag
fdndligo  S190701ah
fdndligo  S190706ai

justt02   GW151226
justt02   GW170104
justt02   GW170608
justt02   GW170729
justt02   GW170809
justt02   GW170814
justt02   GW170817
justt02   GW170818
justt02   GW170823
justt02   S190412m
justt02   S190421ar
justt02   S190425z
justt02   S190503bf
justt02   S190510g
justt02   S190512at
justt02   S190517h
justt02   S190519bj
justt02   S190521r
justt02   S190707q
