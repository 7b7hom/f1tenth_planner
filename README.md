# F1Tenth_planner
Local Planner

## gen_spline
2025.07.07
- calcSpline()까지 구현

2025.07.08
- visual 전까지 완료(orin 보드를 못 사용해서 시각화 전까지의 로직 구현)

2025.07.08
- visual, runPlanningPipeline ai 돌려서 대충 채워놓고, gen_spline.cpp에서 모든 로직 실행하도록 정리함.(실행 목적..)

2025.07.10
- 빌드 및 실행 시도

2025.07.11
- debuging 완료
    -> 실행 잘 됨!
- gui 안 뜨는 문제 해결
  1. VcXcrv Windows X Server install
  2. 방화벽 설정에서 X Server 허용
  3. echo DISPLAY=0.0 두고 실행
 
<img width="640" height="480" alt="graph" src="https://github.com/user-attachments/assets/a0a209e5-5e37-4a7d-b205-5c66fc85abd1" />


## gen_spline 흐름
1. computeEuclideanDistances
  - 주어진 경로(path) 상의 연속된 점들 간의 유클리드 거리 계산
  - 입력: path(MatrixXd 형태의 경로)
  - 출력: 인접한 지점들 사이의 거리 담은 VectorXd

2. calcSplines
  - 주어진 경로에 대한 3차 spline 계수 계산
  - 입력: path(MatrixXd 형태의 경로), el_lengths_ptr(선택적으로 미리 계산된 각 구간의 길이를 가리키는 포인터), psi_s & psi_e: spline의 시작점과 끝점의 heading 제약 조건, use_dist_scaling(미분 연속성 조건에 구간 길이 기반의 스케일링을 적용할지 여부)
  - 출력: SplineResult 객체(계수 행렬 coeffs_x, coeffs_y, 시스템 행렬 M, 정규화된 법선 벡터 normvec_normalized, 구간 길이 ds)
  - 입력된 정보를 기반으로 시스템 행렬 M과 우변 벡터 b_x, b_y 구성(주어진 두 점을 연결하는 3차 스플라인의 수학적 계수 계산)
    1. 연립 선형 방정식(M⋅coeffs=b) 구성: 제약 조건등을 M 행렬과 b_x, b_y 벡터에 채워넣음.
    2. 커리 스케일링: ds 벡터를 기반으로 scaling 팩터 계산하여 미분 연속성 조건에 적용
    3. 방정식 풀이: Eigen::colPivHouseholderQr().solve() 사용하여 coeffs_x, coeffs_y 계산
    4. 법선 벡터 계산: coeffs_x, coeffs_y 1차항 계수(a1, b1) 이용하여 각 spline의 법선 벡터 계산, 정규화
  - generateGraphEdges()에서 각 노드 쌍에 대해 호출됨.
- 
    
