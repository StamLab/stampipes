BEGIN{
  readsep = "~";
}
{
  split($0, a, readsep);
  print a[1];
  print a[2]
}
