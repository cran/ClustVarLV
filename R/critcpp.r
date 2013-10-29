critcpp = function(a,b,c) {
  res = .Call("critcpp",a , b, c,
              package="CLV")
  return(res)
}