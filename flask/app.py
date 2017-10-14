from flask import Flask, jsonify,request, abort
from cproj.tasks import solveproblem
from celery import group

app = Flask(__name__)

PROBLEMS = ["prob1aI", "prob1bI", "prob1cI","prob1bII"]

#PROBLEMS = ["prob1aI", "prob1cI"]

METHODS = ['MC','MC-S','QMC-S','MLMC','MLMC-A',
    'FFT','FGL','COS',
    'FD','FD-NU','FD-AD',
    'RBF','RBF-FD','RBF-PUM','RBF-LSML','RBF-AD','RBF-MLT']

PARAM_NAMES = ['S','K','T','r','sig']

@app.route('/benchop/api/allprobs', methods=['GET'])
def solve_all_problems():
    results = group(solveproblem.s(PROBLEMS[i]) for i in xrange(len(PROBLEMS)))().get()
    ret_dic ={}
    for i in range(len(PROBLEMS)):
        time , err = clean_res(results[i])
        ret_dic[PROBLEMS[i]] = make_ret_dic(time,err,METHODS)
    return jsonify(ret_dic)


@app.route('/benchop/api/prob/<problem_name>', methods=['POST'])
def solve_problem(problem_name):
    if not request.json or not all(par in request.json for par in PARAM_NAMES):
        abort(400)
    result = solveproblem.delay(problem_name,request.json).get()
    time, err = clean_res(result)
    ret = make_ret_dic(time,err,METHODS)
    ret['par'] = request.json
    return jsonify(ret)


def make_rank_dic(m,dic,methods):
    ret_l = [float("inf")]
    p = -1
    while len(ret_l) < len(dic):
        min_v = float("inf")
        for i in range(len(dic)):
            if dic[methods[i]][m] < min_v and dic[methods[i]][m] > min(ret_l):
                p = i
        ret_l.append({methods[p] : methods[p][m]})
    ret_l.pop(0)
    return ret_l



def make_ret_dic(time,err,methods):
    dic={}
    for j in range(len(methods)):
        dic[methods[j]] = {"time" : time[j] , "err" : err[j]}
    return dic


def clean_res(res):
    time_res = []
    err_res = []
    for i in range(len(res[0])):
        time_res.append(res[1][i][0])
        err_res.append(res[0][i][0])
    return time_res, err_res

if __name__ == '__main__':
     app.run(host='0.0.0.0',debug=True)
