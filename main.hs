#!/usr/bin/env nix-shell
#!nix-shell --pure -i runghc -p "haskellPackages.ghcWithPackages (p: [ p.mwc-random p.pointless-fun])"

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
import Data.Function.Pointless

type Point = [Double]
newtype Cluster = UnsafeMkCluster ([Point], [(Point, Double)])
type Clustering = [Cluster]

clustering :: Int -> Double -> [[Point]] -> IO Clustering
clustering t δ v = sequence $ map (cluster (length v) t δ) v

sil :: Clustering -> Double
sil self = let
    my_map :: ([a] -> [a] -> a -> b) -> [a] -> [a] -> [b]
    my_map _ _ [] = []
    my_map f prev (c:succ) = f prev succ c : my_map f (c:prev) succ
  in sum (my_map (sum_sil .: (++)) [] self)
    / fromIntegral (sum $ map (length . points) self)

cluster :: Int -> Int -> Double -> [Point] -> IO Cluster
cluster k t δ c = let
  c_sampled = if length c <= t then return $ zip c [1, 1..] else let
    unif = 1 / fromIntegral (length c)
    prob = 2 * unif * log (2 * fromIntegral k / δ)
    in do
      gen <- createSystemRandom
      let ber_gen p = bernoulli p gen
      s_0 <- filterM (const $ ber_gen prob) c
      let w_0 = map (sum . flip map c . d) s_0
          probs = map (min 1 . (fromIntegral t *) . max unif . maximum . \p -> zipWith ((/) . d p) s_0 w_0) c
      filterM (ber_gen . snd) $ zip c probs
  in UnsafeMkCluster <$> (,) c <$> c_sampled

points :: Cluster -> [Point]
points (UnsafeMkCluster (c, _)) = c

sum_sil :: [Cluster] -> Cluster -> Double
sum_sil others self = sum $ map (sil_point self others) $ points self

w :: Point -> Cluster -> Double
w p (UnsafeMkCluster (_, self)) = sum $ map (uncurry $ (/) . d p) self

sil_point :: Cluster -> [Cluster] -> Point -> Double
sil_point c others self = let
  a = w self c / (fromIntegral . pred . length . points) c
  b = foldl (min .^ liftM2 (/) (w self) (fromIntegral . length . points))
    (1/0) -- infinite
    others
  in (b - a) / max a b

d :: Point -> Point -> Double
d = sqrt . sum .: zipWith ((**2) .: (-))

main :: IO ()
main = let
  c = clustering 4 0.1 [
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
    ]
  in print =<< show <$> sil <$> c

