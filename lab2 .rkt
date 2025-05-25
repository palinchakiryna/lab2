#lang racket

;; --- Bessel I0 approximation ---
(define (bessel-i0 x)
  (define terms 20)
  (define (fact n)
    (if (= n 0) 1 (* n (fact (- n 1)))))
  (define (term k)
    (define num (expt (* 0.5 x) k))
    (/ (sqr num) (sqr (fact k))))
  (for/sum ([k (in-range terms)])
    (term k)))

;; --- Rice PDF ---
(define (rice-pdf x ν σ)
  (define z (/ (* x ν) (* σ σ)))
  (define b (bessel-i0 z))
  (* (/ x (* σ σ))
     (exp (- (/ (+ (* x x) (* ν ν)) (* 2 σ σ))))
     b))

;; --- CDF via trapezoidal rule ---
(define (rice-cdf x ν σ)
  (define N 1000)
  (define dx (/ x N))
  (for/fold ([acc 0.0]) ([i (in-range N)])
    (define xi (* i dx))
    (define xi+1 (* (+ i 1) dx))
    (+ acc (* 0.5 dx (+ (rice-pdf xi ν σ) (rice-pdf xi+1 ν σ))))))

;; --- Approximate inverse CDF (quantile function) ---
(define (rice-quantile p ν σ [tol 0.0001])
  (define max-x (+ ν (* 6 σ)))
  (define max-iters 100)
  (define (bin-search low high iter)
    (if (>= iter max-iters)
        (/ (+ low high) 2)
        (let* ([mid (/ (+ low high) 2)]
               [c (rice-cdf mid ν σ)])
          (cond
            [(not (real? c)) (/ (+ low high) 2)]
            [(< (abs (- c p)) tol) mid]
            [(< c p) (bin-search mid high (add1 iter))]
            [else (bin-search low mid (add1 iter))]))))
  (bin-search 0 max-x 0))

;; --- Read numbers from CSV column "Price"
(define (read-price-series path)
  (define port (open-input-file path))
  (define lines (port->lines port))
  (define data (cdr lines))
  (define values
    (filter-map
     (lambda (line)
       (define columns (string-split line ","))
       (if (>= (length columns) 2)
           (let* ([raw (list-ref columns 1)]
                  [no-quotes (string-replace raw "\"" "")]
                  [no-commas (string-replace no-quotes "," "")]
                  [num (string->number no-commas)])
             num)
           #f))
     data))
  (close-input-port port)
  values)

(define (get-sorted-range numbers)
  (let ([sorted (sort numbers <)])
    (values sorted (first sorted) (last sorted))))

(define (create-rice-intervals min-val max-val alphabet-size)
  (define ν (/ (- max-val min-val) 4))
  (define σ (/ (- max-val min-val) 6))
  (define probs
    (build-list (add1 alphabet-size)
      (lambda (i) (+ 0.05 (* 0.9 (/ i alphabet-size))))))
  (define quantiles (map (lambda (p) (rice-quantile p ν σ)) probs))
  (define min-q (apply min quantiles))
  (define max-q (apply max quantiles))
  (define scaled
    (map (lambda (q)
           (+ min-val (* (/ (- q min-q) (- max-q min-q))
                         (- max-val min-val))))
         quantiles))
  (printf "Межі інтервалів: ~a\n" scaled)
  scaled)

(define (build-alphabet size)
  (if (<= size 26)
      (map (lambda (i) (string (integer->char (+ 65 i)))) (range size))
      (map (lambda (i) (format "A~a" i)) (range size))))

(define (count-intervals numbers intervals alphabet)
  (define counts (make-vector (length alphabet) 0))
  (for ([num numbers])
    (define idx
      (or (for/first ([i (in-range (sub1 (length alphabet)))]
                      [upper (in-list (rest intervals))]
                      #:when (<= num upper))
            i)
          (sub1 (length alphabet))))
    (vector-set! counts idx (add1 (vector-ref counts idx))))
  (displayln "Розподіл чисел по інтервалах:")
  (for ([i (in-range (length alphabet))])
    (displayln (format "Інтервал ~a: ~a чисел"
                       (list-ref alphabet i)
                       (vector-ref counts i)))))

(define (number->letter num intervals alphabet)
  (define size (length alphabet))
  (or (for/first ([i (in-range size)]
                  [upper (in-list (rest intervals))]
                  #:when (<= num upper))
        (list-ref alphabet i))
      (list-ref alphabet (sub1 size))))

(define (numeric->linguistic numbers intervals alphabet)
  (map (lambda (num) (number->letter num intervals alphabet)) numbers))

(define (in-adjacent-pairs seq)
  (if (< (length seq) 2)
      '()
      (for/list ([i (in-range (- (length seq) 1))])
        (list (list-ref seq i) (list-ref seq (add1 i))))))

(define (build-precedence-matrix ling-seq alphabet)
  (define size (length alphabet))
  (define matrix (make-vector (* size size) 0))
  (for ([pair (in-adjacent-pairs ling-seq)])
    (let ([i (index-of alphabet (first pair))]
          [j (index-of alphabet (second pair))])
      (when (and i j)
        (vector-set! matrix (+ (* i size) j)
                     (add1 (vector-ref matrix (+ (* i size) j)))))))
  matrix)

(define (format-matrix matrix size alphabet)
  (define (row->string row-idx)
    (string-join
     (for/list ([j (in-range size)])
       (number->string (vector-ref matrix (+ (* row-idx size) j))))
     " "))
  (string-join
   (cons
    (string-join (cons "" alphabet) " ")
    (for/list ([i (in-range size)])
      (string-join (list (list-ref alphabet i) (row->string i)) " ")))
   "\n"))

(define (get-alphabet-size)
  (printf "Введіть потужність алфавіту (ціле число >= 2): ")
  (define input (string->number (read-line)))
  (if (and (integer? input) (>= input 2))
      input
      (begin
        (printf "Помилка: введіть ціле число >= 2\n")
        (get-alphabet-size))))

(define (solve-task csv-file)
  (define alphabet-size (get-alphabet-size))
  (define alphabet (build-alphabet alphabet-size))
  (define numbers (read-price-series csv-file))
  (when (null? numbers) (error "Список чисел порожній"))
  (define-values (sorted min-val max-val) (get-sorted-range numbers))
  (define intervals (create-rice-intervals min-val max-val alphabet-size))
  (count-intervals numbers intervals alphabet)
  (define ling-seq (numeric->linguistic numbers intervals alphabet))
  (define matrix (build-precedence-matrix ling-seq alphabet))
  (values ling-seq (format-matrix matrix alphabet-size alphabet)))

;; CSV filename
(define csv-file "B-C-D-E-F-Dow Jones Industrial Average Historical Data.csv")

(define-values (ling-seq matrix)
  (let ([init-start-time (current-inexact-milliseconds)])
    (define init-end-time (current-inexact-milliseconds))
    (define exec-start-time (current-inexact-milliseconds))
    (define-values (ling matrix) (solve-task csv-file))
    (define exec-end-time (current-inexact-milliseconds))

    ;; Вывод результатов
    (let ([n 20])
      (printf "Лінгвістичний ряд (перші ~a і останні ~a символів):\n~a ... ~a\n"
              n n
              (take ling (min n (length ling)))
              (take-right ling (min n (length ling)))))
    (printf "Матриця передування:\n~a\n" matrix)
    (printf "Час ініціалізації: ~a сек.\n"
            (/ (- init-end-time init-start-time) 1000.0))
    (printf "Час виконання: ~a сек.\n"
            (/ (- exec-end-time exec-start-time) 1000.0))
    (values ling matrix)))
