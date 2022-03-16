#!/usr/bin/env nix-shell
#!nix-shell --pure -i runghc -p "haskellPackages.ghcWithPackages (p: [ p.mwc-random ])"

--  Copyright 2022 Davide Peressoni
--
--  Licensed under the Apache License, Version 2.0 (the "License");
--  you may not use this file except in compliance with the License.
--  You may obtain a copy of the License at
--
--    http://www.apache.org/licenses/LICENSE-2.0
--
--  Unless required by applicable law or agreed to in writing, software
--  distributed under the License is distributed on an "AS IS" BASIS,
--  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
--  See the License for the specific language governing permissions and
--  limitations under the License.

import System.Random.MWC
import System.Random.MWC.Distributions
import Control.Monad

type Point = [Double]
newtype Cluster = UnsafeMkCluster ([Point], [(Point, Double)])
type Clustering = [Cluster]

clustering :: [[Point]] -> Int -> Double -> IO Clustering
clustering v t δ = let k = length v in
  sequence $ map (\x -> cluster x k t δ) v

sil :: Clustering -> Double
sil self = let
    my_map :: ([a] -> a -> [a] -> b) -> [a] -> [a] -> [b]
    my_map _ _ [] = []
    my_map f prev (c:succ) = f prev c succ : my_map f (c:prev) succ
  in sum (my_map (\prev c succ -> sum_sil c $ prev ++ succ) [] self)
    / fromIntegral (sum $ map (\(UnsafeMkCluster(c,_)) -> length c) self)

cluster :: [Point] -> Int -> Int -> Double -> IO Cluster
cluster c k t δ = let
  c_sampled = if length c <= t then return $ zip c [1, 1..] else let
    unif = 1 / (fromIntegral $ length c)
    prob = 2 * unif * log (2 * (fromIntegral k) / δ)
    in do
      gen <- createSystemRandom
      let ber_gen p = bernoulli p gen
      s_0 <- filterM (\_ -> ber_gen prob) c
      let w_0 = map (\p -> sum $ map (d p) c) s_0
          probs = map (\p -> min 1 (fromIntegral t * max unif (maximum $ map (\(e, w) -> d p e / w) $ zip s_0 w_0)) :: Double) c
      filterM (ber_gen . snd) $ zip c probs
  in UnsafeMkCluster <$> (,) c <$> c_sampled

sum_sil :: Cluster -> [Cluster] -> Double
sum_sil self@(UnsafeMkCluster (c, _)) others = sum $ map (\e -> sil_point e self others) c

w :: Cluster -> Point -> Double
w (UnsafeMkCluster (_, self)) p = sum $ map (\(e, p_e) -> d p e / p_e) self

sil_point :: Point -> Cluster -> [Cluster] -> Double
sil_point self this_c others = let
  a = w this_c self / fromIntegral (let UnsafeMkCluster(c, _) = this_c in length c - 1)
  b = foldl (\m c -> min m $ w c self / fromIntegral (let UnsafeMkCluster(o, _) = c in length o))
    (1/0) -- infinite
    others
  in (b - a) / max a b

d :: Point -> Point -> Double
d self p = sqrt $ sum $ (**2) <$> zipWith (-) self p

main :: IO ()
main = let
  c = clustering [
        -- exact silhouette 0.19340581656247494
                   [
                     [1, 1],
                     [2, 2],
                     [3, 2],
                     [4, 2],
                     [-1, 3],
                     [5, -1],
                     [60.3, 24.5],
                     [1, 0.5],
                     [5, 3],
                     [20, 21]
                   ],
                   [
                     [-1, -1],
                     [-2, -2],
                     [-3, -2],
                     [-4, -2],
                     [1, -3],
                     [-5, 1],
                     [-60.3, -24.5],
                     [-1, -0.5],
                     [-5, -3],
                     [-20, -21]
                   ]
                 ] 4 0.1
  in show <$> sil <$> c >>= print

