from flask import Flask, jsonify
from cproj.tasks import solveproblem
from celery import group

app = Flask(__name__)

PROBLEMS = ["prob1aI", "prob1bI", "prob1cI", "prob1aII","prob1bII", "prob1cII"]

METHODS = ['MC','MC-S','QMC-S','MLMC','MLMC-A',
    'FFT','FGL','COS',
    'FD','FD-NU','FD-AD',
    'RBF','RBF-FD','RBF-PUM','RBF-LSML','RBF-AD','RBF-MLT']



@app.route('/benchop/api/allprobs', methods=['GET'])
def solve_all_problems():
    results = group(solveproblem.s(PROBLEMS[i]) for i in xrange(len(PROBLEMS)))().get()
    ret_dic ={}
    for i in range(len(PROBLEMS)):
        time , err = clean_res(result[i])
        dic ={}
        for j in range(len(METHODS)):
            dic[METHODS[j]] = {"time" : time[j] , "err" : err[j]}
        ret_dic[PROBLEMS[i]] = dic
    return jsonify(ret_dic)

def clean_res(res):
    time_res = []
    err_res = []
    for i in range(len(res[0])):
        time_res.append(res[0][i][0])
        err_res.append(res[1][i][0])
    return time_res, err_res

if __name__ == '__main__':
     app.run(host='0.0.0.0',debug=True)
