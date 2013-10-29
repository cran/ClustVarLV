mincpp = function(a) {
  res = .Call("mincpp",a ,
              package="CLV")
  return(res)
}