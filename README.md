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

## gen_spline 흐름
- computeEuclideanDistances
  - 주어진 경로(path) 상의 연속된 점들 간의 유클리드 거리 계산
  - 입력: path(MatrixXd 형태의 경로) / 출력: 인접한 지점들 사이의 거리 담은 VectorXd

- calcSplines(핵심)
  - 주어진 경로에 대한 3차 spline 계수 계산
  - 입력: path(MatrixXd 형태의 경로), el_lengths_ptr(선택적으로 미리 계산된 각 구간의 길이를 가리키는 포인터), psi_s & psi_e: spline의 시작점과 끝점의 heading 제약 조건, use_dist_scaling(미분 연속성 조건에 구간 길이 기반의 스케일링을 적용할지 여부)
  - 입력된 정보를 기반으로 시스템 행렬 M과 우변 벡터 b_x, b_y 구성  
